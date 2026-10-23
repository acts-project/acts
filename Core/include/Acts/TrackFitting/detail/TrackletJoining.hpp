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
#include "Acts/EventData/TrackProxyCommon.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/TrackStateProxyConcept.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cstddef>
#include <optional>
#include <system_error>
#include <utility>
#include <vector>

namespace Acts::detail {

/// One smoothing step from the first state of a smoothed segment onto the
/// prior node it was predicted from, in the Rauch-Tung-Striebel gain form of
/// the @ref Acts::GainMatrixSmoother: `G = P F^T Pp^-1`,
/// `x' = x + G (xs - xp)`, `P' = P + G (Ps - Pp) G^T`
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
///         covariance is not positive definite
Result<SharedStateCondition> makeSharedStateCondition(
    const AnyConstTrackStateProxy& shared, const BoundVector& newMean,
    const BoundMatrix& newCov, const Logger& logger);

/// Sensitivity of the filtered state to the prior mean, `Df = (I - K H) Dp`.
/// This and @ref conditionTrackState are the recursions of the modified
/// Bryson-Frazier smoother (@ref Acts::MbfSmoother) applied to the
/// sensitivity to the prior mean instead of to the innovations.
/// @param ts The state
/// @param dp The sensitivity of its predicted state
/// @return The sensitivity of its filtered state
BoundMatrix filteredSensitivity(const AnyConstTrackStateProxy& ts,
                                const BoundMatrix& dp);

/// Condition one smoothed state and transport the adjoint sensitivity to the
/// previous state, as the modified Bryson-Frazier smoother transports its
/// adjoint: `T = Df - Pf Psi^`, `x += T u`, `P += T w T^T`,
/// `Psi^ <- F^T (H^T S^-1 H Dp + (I - K H)^T Psi^)`
/// @param ts The state
/// @param dp The sensitivity of its predicted state
/// @param df The sensitivity of its filtered state
/// @param psiHat The adjoint sensitivity, in and out
/// @param condition The per-track terms
/// @param logger The logger
/// @return Success, or `KalmanFitterError::SmoothFailed` if the result is
///         not finite
Result<void> conditionTrackState(AnyMutableTrackStateProxy ts,
                                 const BoundMatrix& dp, const BoundMatrix& df,
                                 BoundMatrix& psiHat,
                                 const SharedStateCondition& condition,
                                 const Logger& logger);

/// Condition the smoothed states of a segment on a new estimate of one of
/// its states.
///
/// The segment runs from @p sharedIndex to @p tipIndex in filter order and is
/// smoothed. Each state is updated by its regression on the shared state,
/// `x_k += L_k (x_s' - x_s)`, `P_k += L_k (P_s' - P_s) L_k^T`, which is the
/// exact linear-Gaussian update for information that enters through the
/// shared state. `L_k` is the sensitivity of the smoothed state to the prior
/// mean of the segment, computed in the modified Bryson-Frazier form with a
/// forward sweep over @ref filteredSensitivity and a backward sweep over
/// @ref conditionTrackState.
///
/// States before the shared state get a plain smoothing step from the new
/// estimate. States after the last smoothed state stay as they are.
///
/// @tparam traj_t The MultiTrajectory type
/// @param trajectory The trajectory that holds the segment
/// @param sharedIndex The index of the shared state
/// @param tipIndex The index of the last state of the segment in filter order
/// @param newMean The new estimate of the shared state
/// @param newCov Its covariance
/// @param logger The logger
/// @return Success, or `KalmanFitterError::SmoothFailed` if the shared
///         state is not on the chain, a state lacks components, or a result
///         is not finite
template <typename traj_t>
Result<void> conditionSmoothedStatesOnSharedState(
    traj_t& trajectory, std::size_t sharedIndex, std::size_t tipIndex,
    const BoundVector& newMean, const BoundMatrix& newCov,
    const Logger& logger = getDummyLogger()) {
  using TrackStateProxy = typename traj_t::TrackStateProxy;

  // collect the segment from the tip down to the shared state
  std::vector<std::size_t> indices;
  bool found = false;
  trajectory.applyBackwards(tipIndex, [&](TrackStateProxy ts) {
    indices.push_back(ts.index());
    found = (ts.index() == sharedIndex);
    return !found;
  });
  if (!found) {
    ACTS_DEBUG("Shared state " << sharedIndex << " is not on the chain of "
                               << tipIndex);
    return KalmanFitterError::SmoothFailed;
  }

  // states after the last smoothed state, e.g. holes after the last
  // measurement, stay as they are
  std::size_t skip = 0;
  while (skip + 1 < indices.size() &&
         !trajectory.getTrackState(indices[skip]).hasSmoothed()) {
    ++skip;
  }
  indices.erase(indices.begin(), indices.begin() + skip);

  auto shared = trajectory.getTrackState(sharedIndex);
  if (!shared.hasSmoothed()) {
    ACTS_DEBUG("Shared state " << sharedIndex << " has no smoothed parameters");
    return KalmanFitterError::SmoothFailed;
  }
  auto condition = makeSharedStateCondition(AnyConstTrackStateProxy{shared},
                                            newMean, newCov, logger);
  if (!condition.ok()) {
    return condition.error();
  }

  // forward sweep: sensitivities (Dp, Df) of the predicted and filtered
  // states to the prior mean of the segment
  std::vector<std::pair<BoundMatrix, BoundMatrix>> sensitivities(
      indices.size());
  BoundMatrix df = BoundMatrix::Identity();
  for (std::size_t k = indices.size(); k-- > 0;) {
    auto ts = trajectory.getTrackState(indices[k]);
    if (!ts.hasSmoothed() || !ts.hasPredicted() || !ts.hasJacobian()) {
      ACTS_DEBUG("State " << ts.index()
                          << " lacks smoothed, predicted or jacobian");
      return KalmanFitterError::SmoothFailed;
    }
    const bool isShared = (k + 1 == indices.size());
    const BoundMatrix dp =
        isShared ? BoundMatrix::Identity() : BoundMatrix(ts.jacobian() * df);
    df = filteredSensitivity(AnyConstTrackStateProxy{ts}, dp);
    sensitivities[k] = {dp, df};
  }

  // backward sweep: condition the states from the tip down to the shared one
  BoundMatrix psiHat = BoundMatrix::Zero();
  for (std::size_t k = 0; k + 1 < indices.size(); ++k) {
    auto ts = trajectory.getTrackState(indices[k]);
    const auto& [dp, dfk] = sensitivities[k];
    auto res = conditionTrackState(AnyMutableTrackStateProxy{ts}, dp, dfk,
                                   psiHat, *condition, logger);
    if (!res.ok()) {
      return res.error();
    }
  }

  shared.smoothed() = newMean;
  shared.smoothedCovariance() = newCov;

  if (shared.hasPrevious()) {
    GainMatrixSmoother smoother;
    std::error_code error;
    auto later = shared;
    trajectory.applyBackwards(shared.previous(), [&](TrackStateProxy ts) {
      ts.addComponents(TrackStatePropMask::Smoothed);
      auto res = smoother.calculate(AnyMutableTrackStateProxy{ts},
                                    AnyConstTrackStateProxy{later}, logger);
      if (!res.ok()) {
        error = res.error();
        return false;
      }
      later = ts;
      return true;
    });
    if (error) {
      return Result<void>::failure(error);
    }
  }

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
                                              bool invertJacobians = true) {
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
