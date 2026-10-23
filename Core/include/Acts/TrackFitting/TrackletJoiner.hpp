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
#include "Acts/TrackFitting/SmootherAdjoint.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <cstddef>
#include <optional>
#include <tuple>

namespace Acts {

/// @addtogroup track_fitting
/// @{

/// Join a second tracklet to a first one through a shared state, and give
/// every state the information of both.
///
/// The first tracklet is smoothed and holds its own measurements. The second
/// tracklet was filtered from one state of the first, the shared state, as
/// prior; its filtered states hold the first tracklet plus its own
/// measurements up to each state. This is the two-way track finding, in either
/// direction, but the operation is general.
///
/// A smoother over the stitched chain would count the first tracklet twice.
/// The joiner instead
///
/// 1. smooths the second tracklet alone and carries the adjoint onto the
///    shared state (@ref SmootherAdjoint), which gives the shared state the
///    information of both tracklets,
/// 2. transports that change to every state of the first tracklet by
///    conditioning (@ref conditionSmoothedStatesOnSharedState), which is exact
///    in the linear-Gaussian model and needs no propagation,
/// 3. links the second tracklet behind the shared state
///    (@ref linkTrackletBehind).
///
/// Afterwards every state holds smoothed parameters with the information of
/// both tracklets. The predicted and filtered history stays in two filter
/// orders, so a later smoother over the joined chain is invalid. The track
/// level data of the joined track is taken from the second tracklet, which
/// holds the counts of both and, where the second pass extrapolated, the
/// parameters at the reference surface.
struct TrackletJoiner {
  /// Configuration
  struct Config {
    /// Smooth the second tracklet and carry its information onto the shared
    /// state
    bool smoothSecond = true;
    /// Transport the change of the shared state to the first tracklet
    bool conditionFirst = true;
    /// Invert and shift the jacobians of the second tracklet in the link
    bool invertJacobians = true;
  };

  /// The configuration
  Config config;

  /// Default constructor
  TrackletJoiner() = default;

  /// Constructor
  /// @param cfg The configuration
  explicit TrackletJoiner(const Config& cfg) : config(cfg) {}

  /// Join in place. The first tracklet becomes the joined track; the second
  /// tracklet's states move into its chain.
  ///
  /// The first tracklet must own its states: with several second tracklets
  /// from one first tracklet, join copies with @ref joinCopy.
  ///
  /// @tparam smoother_t The smoother for the second tracklet
  /// @tparam track_proxy_t The track proxy type
  /// @param gctx The geometry context
  /// @param first The first tracklet, smoothed
  /// @param sharedState The state of the first tracklet the second started from
  /// @param second The second tracklet, filtered from @p sharedState
  /// @param logger The logger
  /// @param smoother The smoother
  /// @return Success or the first error
  template <typename smoother_t = MbfSmoother, TrackProxyConcept track_proxy_t>
  Result<void> joinInPlace(const GeometryContext& gctx, track_proxy_t& first,
                           typename track_proxy_t::TrackStateProxy sharedState,
                           track_proxy_t& second,
                           const Logger& logger = getDummyLogger(),
                           const smoother_t& smoother = smoother_t{}) const {
    auto& trajectory = first.container().trackStateContainer();

    if (second.tipIndex() == kTrackIndexInvalid) {
      ACTS_DEBUG("Second tracklet has no states");
      return KalmanFitterError::SmoothFailed;
    }
    if (!sharedState.hasSmoothed()) {
      ACTS_DEBUG("Shared state " << sharedState.index()
                                 << " has no smoothed parameters");
      return KalmanFitterError::SmoothFailed;
    }

    BoundVector sharedMean = sharedState.smoothed();
    BoundMatrix sharedCov = sharedState.smoothedCovariance();

    if (config.smoothSecond) {
      SmootherAdjoint adjoint;
      auto res = smoother.smooth(gctx, trajectory, second.tipIndex(),
                                 kTrackIndexInvalid, &adjoint, logger);
      if (!res.ok()) {
        ACTS_DEBUG("Smoothing of the second tracklet failed: " << res.error());
        return res.error();
      }
      std::tie(sharedMean, sharedCov) =
          applySmootherAdjoint(sharedMean, sharedCov, adjoint);
    }

    if (config.conditionFirst) {
      auto res = conditionSmoothedStatesOnSharedState(
          trajectory, sharedState.index(), first.tipIndex(), sharedMean,
          sharedCov, logger);
      if (!res.ok()) {
        ACTS_DEBUG(
            "Conditioning of the first tracklet failed: " << res.error());
        return res.error();
      }
    } else {
      sharedState.smoothed() = sharedMean;
      sharedState.smoothedCovariance() = sharedCov;
    }

    const auto tipIndex = first.tipIndex();
    linkTrackletBehind(sharedState, second, config.invertJacobians);
    first.copyFromWithoutStates(second);
    first.tipIndex() = tipIndex;
    first.stemIndex() = second.stemIndex();

    return Result<void>::success();
  }

  /// Join copies. Both tracklets are copied into @p out, the copies are
  /// joined, and the joined copy is returned. The inputs stay unchanged, so
  /// several second tracklets can be joined to one first tracklet.
  ///
  /// @tparam smoother_t The smoother for the second tracklet
  /// @tparam track_container_t The output track container type
  /// @tparam track_proxy_t The input track proxy type
  /// @tparam track_state_proxy_t The input track state proxy type
  /// @param gctx The geometry context
  /// @param out The container that receives the copies
  /// @param first The first tracklet, smoothed
  /// @param sharedState The state of the first tracklet the second started from
  /// @param second The second tracklet, filtered from @p sharedState
  /// @param logger The logger
  /// @param smoother The smoother
  /// @return The joined track in @p out, or the first error
  template <typename smoother_t = MbfSmoother,
            TrackContainerFrontend track_container_t,
            TrackProxyConcept track_proxy_t,
            TrackStateProxyConcept track_state_proxy_t>
  Result<typename track_container_t::TrackProxy> joinCopy(
      const GeometryContext& gctx, track_container_t& out,
      const track_proxy_t& first, const track_state_proxy_t& sharedState,
      const track_proxy_t& second, const Logger& logger = getDummyLogger(),
      const smoother_t& smoother = smoother_t{}) const {
    using TrackProxy = typename track_container_t::TrackProxy;
    using TrackStateProxy = typename track_container_t::TrackStateProxy;

    // position of the shared state, counted from the tip
    std::size_t steps = 0;
    bool found = false;
    for (const auto& ts : first.trackStatesReversed()) {
      if (ts.index() == sharedState.index()) {
        found = true;
        break;
      }
      ++steps;
    }
    if (!found) {
      ACTS_DEBUG("Shared state " << sharedState.index()
                                 << " is not on the first tracklet");
      return Result<TrackProxy>::failure(KalmanFitterError::SmoothFailed);
    }

    TrackProxy firstCopy = out.makeTrack();
    firstCopy.copyFrom(first);
    TrackProxy secondCopy = out.makeTrack();
    secondCopy.copyFrom(second);

    std::optional<TrackStateProxy> sharedCopy;
    std::size_t i = 0;
    for (auto ts : firstCopy.trackStatesReversed()) {
      if (i == steps) {
        sharedCopy = ts;
        break;
      }
      ++i;
    }

    auto res =
        joinInPlace(gctx, firstCopy, *sharedCopy, secondCopy, logger, smoother);
    // the second copy is a shell now: its states live in the joined chain
    out.removeTrack(secondCopy.index());
    if (!res.ok()) {
      return Result<TrackProxy>::failure(res.error());
    }
    return firstCopy;
  }
};

/// @}

}  // namespace Acts
