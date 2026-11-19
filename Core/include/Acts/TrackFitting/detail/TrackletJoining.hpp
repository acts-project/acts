// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/EventData/TrackProxyCommon.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/TrackStateProxyConcept.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/TrackFitting/detail/SmoothedConditioning.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include <Eigen/Cholesky>

namespace Acts::detail {

/// One smoothing step from the first state of a smoothed segment onto the
/// prior node it was predicted from, in the Rauch-Tung-Striebel gain form of
/// the @ref Acts::GainMatrixSmoother: `G = P F^T Pp^-1`,
/// `x' = x + G (xs - xp)`, `P' = P + G (Ps - Pp) G^T`. The gain form needs
/// one factorization of `Pp` per segment. The modified Bryson-Frazier form
/// would need the adjoint of the first state, which the smoother does not
/// store and which costs the same factorization to recover.
/// @param first The first state of the segment, smoothed
/// @param priorMean The estimate at the prior node
/// @param priorCov Its covariance
/// @param logger The logger
/// @return The estimate at the prior node with the segment included, or
///         `KalmanFitterError::SmoothFailed` if the predicted covariance
///         is not positive definite or the result is not finite
Result<std::pair<BoundVector, BoundMatrix>> smoothPriorNode(
    const AnyConstTrackStateProxy& first, const BoundVector& priorMean,
    const BoundMatrix& priorCov, const Logger& logger);

/// Per-track terms of the conditioning on a shared state: the change of the
/// shared state mapped by `Z = Pp_s Ps_s^-1`
struct SharedStateCondition {
  /// `Z (x_s' - x_s)`
  BoundVector u = BoundVector::Zero();
  /// `Z (P_s' - P_s) Z^T`
  BoundMatrix w = BoundMatrix::Zero();
};

/// Build the conditioning terms from the shared state and its new estimate
/// @param shared The shared state, smoothed
/// @param newMean The new estimate of the shared state
/// @param newCov Its covariance
/// @param logger The logger
/// @return The terms, or `KalmanFitterError::SmoothFailed` if the shared
///         covariance is not positive definite or the terms are not finite
Result<SharedStateCondition> makeSharedStateCondition(
    const AnyConstTrackStateProxy& shared, const BoundVector& newMean,
    const BoundMatrix& newCov, const Logger& logger);

/// Condition the smoothed states of a chain on a new estimate of its shared
/// state, the exact linear-Gaussian update for information that enters
/// through that state: @ref conditionStatesAfterNode with `u = Z (x_s' - x_s)`
/// and `W = Z (P_s' - P_s) Z^T`, `Z = Pp_s Ps_s^-1`, and
/// @ref conditionStatesBeforeNode. States after the last smoothed state stay
/// as they are.
///
/// @tparam traj_t The MultiTrajectory type
/// @param trajectory The trajectory that holds the segment
/// @param sharedIndex The index of the shared state
/// @param tipIndex The index of the last state of the chain in filter order
/// @param newMean The new estimate of the shared state
/// @param newCov Its covariance
/// @param logger The logger
/// @return Success, or `KalmanFitterError::SmoothFailed` if the shared
///         state is not on the chain, a state lacks components, or a result
///         is not finite
template <typename traj_t>
Result<void> conditionSmoothedStatesOnSharedState(traj_t& trajectory,
                                                  std::size_t sharedIndex,
                                                  std::size_t tipIndex,
                                                  const BoundVector& newMean,
                                                  const BoundMatrix& newCov,
                                                  const Logger& logger) {
  std::vector<AnyMutableTrackStateProxy> states =
      statesInFilterOrder(trajectory, tipIndex);
  const auto it = std::ranges::find_if(
      states, [&](const auto& ts) { return ts.index() == sharedIndex; });
  if (it == states.end()) {
    ACTS_DEBUG("Shared state " << sharedIndex << " is not on the chain of "
                               << tipIndex);
    return KalmanFitterError::SmoothFailed;
  }
  const auto pos = static_cast<std::size_t>(std::distance(states.begin(), it));
  AnyMutableTrackStateProxy shared = *it;
  if (!shared.hasSmoothed() || !shared.hasPredicted()) {
    ACTS_DEBUG("Shared state " << sharedIndex
                               << " has no smoothed or predicted parameters");
    return KalmanFitterError::SmoothFailed;
  }

  auto sharedProxy = trajectory.getTrackState(sharedIndex);
  auto condition = makeSharedStateCondition(
      AnyConstTrackStateProxy{sharedProxy}, newMean, newCov, logger);
  if (!condition.ok()) {
    return condition.error();
  }
  auto res = conditionStatesAfterNode<eBoundSize>(
      std::span(states).subspan(pos + 1),
      filteredSensitivity<eBoundSize>(shared, BoundMatrix::Identity()),
      condition->u, condition->w, logger);
  if (!res.ok()) {
    return res.error();
  }

  if (pos > 0) {
    const Eigen::LDLT<BoundMatrix> ldlt(shared.predictedCovariance());
    if (ldlt.info() != Eigen::Success || !(ldlt.vectorD().array() > 0).all()) {
      ACTS_DEBUG(
          "Predicted covariance of the shared state is not positive "
          "definite");
      return KalmanFitterError::SmoothFailed;
    }
    res = conditionStatesBeforeNode<eBoundSize>(
        std::span(states).first(pos + 1), ldlt.solve(BoundMatrix::Identity()),
        subtractBoundParameters(newMean, shared.smoothed()),
        newCov - shared.smoothedCovariance(), logger);
    if (!res.ok()) {
      return res.error();
    }
  }

  setSmoothedWithoutSharing(shared, newMean, newCov);
  return Result<void>::success();
}

/// Link a second tracklet, filtered from a state of a first tracklet, behind
/// that state.
///
/// The states of the second tracklet are reversed and linked behind the shared
/// state. With @p invertJacobians the jacobians are inverted and shifted, and
/// the transport onto the shared state is stored on the shared state.
///
/// @tparam track_state_proxy_t The track state proxy type
/// @tparam track_proxy_t The track proxy type
/// @param sharedState The shared state of the first tracklet
/// @param second The second tracklet
/// @param invertJacobians Whether to invert and shift the jacobians
/// @return The transport onto the shared state, if any
template <TrackStateProxyConcept track_state_proxy_t,
          TrackProxyConcept track_proxy_t>
std::optional<BoundMatrix> linkTrackletBehind(track_state_proxy_t sharedState,
                                              track_proxy_t& second,
                                              bool invertJacobians) {
  std::optional<BoundMatrix> stitch =
      second.reverseTrackStates(invertJacobians);
  if (second.tipIndex() == kTrackIndexInvalid) {
    return stitch;
  }

  auto secondTip = second.outermostTrackState();
  sharedState.previous() = secondTip.index();
  secondTip.template component<TrackIndexType>(detail_tp::kNextKey) =
      sharedState.index();

  if (stitch.has_value()) {
    sharedState.jacobian() = *stitch;
  }

  return stitch;
}

}  // namespace Acts::detail
