// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackContainerFrontendConcept.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStateProxyConcept.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/MbfSmoother.hpp"
#include "Acts/TrackFitting/detail/TrackletJoining.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <iterator>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace Acts {

/// @addtogroup track_fitting
/// @{

/// Join a second tracklet, filtered from a state of a smoothed first tracklet
/// as prior, behind that shared state, and combine their smoothed estimates
/// without a new propagation. This is the two-way track finding.
///
/// A smoother over the stitched chain would count the first tracklet twice.
/// Instead the second tracklet is smoothed alone and one smoothing step
/// carries its information onto the shared state; the change of the shared
/// state is transported to the first tracklet by conditioning
/// (`detail::conditionSmoothedStatesOnSharedState`); then the chains are linked
/// (`detail::linkTrackletBehind`).
///
/// With smoothing and conditioning enabled, the smoothed states hold both
/// tracklets. Unsmoothed trailing states of the first tracklet stay unchanged.
/// States of the first tracklet before the shared state leave the chain.
/// The predicted and filtered history stays in two filter orders, so a
/// later smoother over the joined chain is invalid. The track-level data is
/// taken from the second tracklet.
/// If the second tracklet has no states, only its track-level data is copied;
/// this covers propagation to the target without another detector surface.
class TrackletJoiner {
 public:
  /// Configuration
  struct Config {
    /// Smooth the second tracklet and update the shared state
    bool smoothSecond = true;
    /// Condition the first tracklet on the updated shared state
    bool conditionFirst = true;
    /// Invert and shift the jacobians of the second tracklet in the link
    bool invertJacobians = true;
  };

  /// The prior a second tracklet was filtered from when it is not the
  /// smoothed shared state, e.g. that state with an inflated covariance. The
  /// join subtracts it at the shared state and conditions the second tracklet
  /// on the result. Without `smoothSecond` it has no effect.
  struct SecondPrior {
    /// The prior parameters
    BoundVector mean = BoundVector::Zero();
    /// The prior covariance
    BoundMatrix covariance = BoundMatrix::Identity();
  };

  /// Default constructor
  TrackletJoiner() : TrackletJoiner(Config{}) {}

  /// Constructor
  /// @param cfg The configuration
  /// @param logger The logger
  explicit TrackletJoiner(const Config& cfg,
                          std::unique_ptr<const Logger> logger =
                              getDefaultLogger("TrackletJoiner", Logging::INFO))
      : m_cfg(cfg), m_logger(std::move(logger)) {}

  /// Get the configuration
  /// @return The configuration
  const Config& config() const { return m_cfg; }

  /// Join in place: the first tracklet becomes the joined track. It must own
  /// its states; with several second tracklets use @ref joinCopy.
  /// Both tracklets must use the same trajectory and own disjoint states.
  /// If the first tracklet is not linked forward, this links it forward.
  /// On failure, smoothed states may already have been modified.
  ///
  /// @tparam smoother_t The smoother for the second tracklet
  /// @tparam track_proxy_t The track proxy type
  /// @param gctx The geometry context
  /// @param first The first tracklet, smoothed
  /// @param sharedState The state the second tracklet started from
  /// @param second The second tracklet
  /// @param smoother The smoother
  /// @param secondPrior The prior of the second tracklet, if not the
  ///        smoothed shared state
  /// @return Success or the first error
  template <typename smoother_t = MbfSmoother, TrackProxyConcept track_proxy_t>
  Result<void> joinInPlace(
      const GeometryContext& gctx, track_proxy_t& first,
      typename track_proxy_t::TrackStateProxy sharedState,
      track_proxy_t& second, const smoother_t& smoother = smoother_t{},
      const std::optional<SecondPrior>& secondPrior = std::nullopt) const {
    auto& trajectory = first.container().trackStateContainer();

    if (first.stemIndex() == kTrackIndexInvalid) {
      first.linkForward();
    }

    if (second.tipIndex() == kTrackIndexInvalid) {
      // An empty second pass can still carry the estimate at the target.
      const auto tipIndex = first.tipIndex();
      const auto stemIndex = first.stemIndex();
      first.copyFromWithoutStates(second);
      first.tipIndex() = tipIndex;
      first.stemIndex() = stemIndex;
      return Result<void>::success();
    }
    if (!sharedState.hasSmoothed()) {
      ACTS_DEBUG("Shared state " << sharedState.index()
                                 << " has no smoothed parameters");
      return KalmanFitterError::SmoothFailed;
    }

    BoundVector sharedMean = sharedState.smoothed();
    BoundMatrix sharedCov = sharedState.smoothedCovariance();

    if (m_cfg.smoothSecond) {
      auto res = smoother(gctx, trajectory, second.tipIndex(), logger());
      if (!res.ok()) {
        ACTS_DEBUG("Smoothing of the second tracklet failed: "
                   << res.error().message());
        return res.error();
      }
      // the first state of the second tracklet was predicted from its prior
      auto secondStart = second.outermostTrackState();
      for (auto ts : second.trackStatesReversed()) {
        secondStart = ts;
      }
      const SecondPrior prior =
          secondPrior.value_or(SecondPrior{sharedMean, sharedCov});
      auto node =
          detail::smoothPriorNode(AnyConstTrackStateProxy{secondStart},
                                  prior.mean, prior.covariance, logger());
      if (!node.ok()) {
        return node.error();
      }
      if (!secondPrior.has_value() ||
          (prior.mean == sharedMean && prior.covariance == sharedCov)) {
        // The matching prior cancels exactly; no information subtraction or
        // conditioning of the second tracklet is needed.
        std::tie(sharedMean, sharedCov) = *node;
      } else {
        auto joint = detail::subtractSecondPrior(
            sharedMean, sharedCov, prior.mean, prior.covariance, node->first,
            node->second, logger());
        if (!joint.ok()) {
          return joint.error();
        }
        std::vector<AnyMutableTrackStateProxy> secondStates =
            detail::statesInFilterOrder(trajectory, second.tipIndex());
        auto res2 = detail::conditionSmoothedStatesOnPriorNode(
            secondStates, prior.covariance, node->first, node->second,
            joint->first, joint->second, logger());
        if (!res2.ok()) {
          ACTS_DEBUG("Conditioning of the second tracklet failed: "
                     << res2.error().message());
          return res2.error();
        }
        std::tie(sharedMean, sharedCov) = *joint;
      }
    }

    if (m_cfg.conditionFirst) {
      auto res = detail::conditionSmoothedStatesOnSharedState(
          trajectory, sharedState.index(), first.tipIndex(), sharedMean,
          sharedCov, logger());
      if (!res.ok()) {
        ACTS_DEBUG("Conditioning of the first tracklet failed: "
                   << res.error().message());
        return res.error();
      }
    } else {
      detail::setSmoothedWithoutSharing(AnyMutableTrackStateProxy{sharedState},
                                        sharedMean, sharedCov);
    }

    const auto tipIndex = first.tipIndex();
    detail::linkTrackletBehind(sharedState, second, m_cfg.invertJacobians);
    first.copyFromWithoutStates(second);
    first.tipIndex() = tipIndex;
    first.stemIndex() = second.stemIndex();

    return Result<void>::success();
  }

  /// Join copies of both tracklets in @p out and return the joined copy. The
  /// inputs stay unchanged.
  ///
  /// @tparam smoother_t The smoother for the second tracklet
  /// @tparam track_container_t The output track container type
  /// @tparam track_proxy_t The input track proxy type
  /// @tparam track_state_proxy_t The input track state proxy type
  /// @param gctx The geometry context
  /// @param out The container that receives the copies
  /// @param first The first tracklet, smoothed
  /// @param sharedState The state the second tracklet started from
  /// @param second The second tracklet
  /// @param smoother The smoother
  /// @param secondPrior The prior of the second tracklet, if not the
  ///        smoothed shared state
  /// @return The joined track in @p out, or the first error
  template <typename smoother_t = MbfSmoother,
            TrackContainerFrontend track_container_t,
            TrackProxyConcept track_proxy_t,
            TrackStateProxyConcept track_state_proxy_t>
  Result<typename track_container_t::TrackProxy> joinCopy(
      const GeometryContext& gctx, track_container_t& out,
      const track_proxy_t& first, const track_state_proxy_t& sharedState,
      const track_proxy_t& second, const smoother_t& smoother = smoother_t{},
      const std::optional<SecondPrior>& secondPrior = std::nullopt) const {
    using TrackProxy = typename track_container_t::TrackProxy;
    using TrackStateProxy = typename track_container_t::TrackStateProxy;

    // position of the shared state, counted from the tip
    auto firstStates = first.trackStatesReversed();
    auto sharedIt = std::ranges::find_if(firstStates, [&](const auto& ts) {
      return ts.index() == sharedState.index();
    });
    if (sharedIt == firstStates.end()) {
      ACTS_DEBUG("Shared state " << sharedState.index()
                                 << " is not on the first tracklet");
      return Result<TrackProxy>::failure(KalmanFitterError::SmoothFailed);
    }

    TrackProxy firstCopy = out.makeTrack();
    firstCopy.copyFrom(first);
    TrackProxy secondCopy = out.makeTrack();
    secondCopy.copyFrom(second);

    auto copyStates = firstCopy.trackStatesReversed();
    TrackStateProxy sharedCopy = *std::ranges::next(
        copyStates.begin(),
        std::ranges::distance(firstStates.begin(), sharedIt));

    auto res = joinInPlace(gctx, firstCopy, sharedCopy, secondCopy, smoother,
                           secondPrior);
    // the second copy is a shell now, its states are in the joined chain
    out.removeTrack(secondCopy.index());
    if (!res.ok()) {
      out.removeTrack(firstCopy.index());
      return Result<TrackProxy>::failure(res.error());
    }
    return firstCopy;
  }

 private:
  const Logger& logger() const { return *m_logger; }

  Config m_cfg;
  std::unique_ptr<const Logger> m_logger;
};

/// @}

}  // namespace Acts
