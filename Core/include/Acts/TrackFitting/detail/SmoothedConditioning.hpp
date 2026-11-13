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
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <cstddef>
#include <span>
#include <vector>

namespace Acts::detail {

/// The states of the chain that ends at @p tipIndex, in filter order
/// @param trajectory The trajectory
/// @param tipIndex The last state of the chain
/// @return The states
template <typename traj_t>
std::vector<AnyMutableTrackStateProxy> statesInFilterOrder(
    traj_t& trajectory, std::size_t tipIndex) {
  std::vector<AnyMutableTrackStateProxy> states;
  trajectory.applyBackwards(tipIndex,
                            [&](auto ts) { states.emplace_back(ts); });
  std::ranges::reverse(states);
  return states;
}

/// Write a smoothed estimate without changing aliased filter history.
/// @param ts The state to update
/// @param mean The smoothed parameters
/// @param cov The smoothed covariance
void setSmoothedWithoutSharing(AnyMutableTrackStateProxy ts, BoundVector mean,
                               BoundMatrix cov);

/// Sensitivity of the filtered state to a variable, `Df = (I - K H) Dp`.
/// @param ts The state
/// @param dp The sensitivity of its predicted state
/// @return The sensitivity of its filtered state
template <std::size_t C>
Matrix<eBoundSize, C> filteredSensitivity(const AnyMutableTrackStateProxy& ts,
                                          const Matrix<eBoundSize, C>& dp);

/// Condition the smoothed states after a node on a change `(u, W)` of a
/// variable: `x_j += T_j u`, `P_j += T_j W T_j^T`, with `T_j` the sensitivity
/// of the smoothed state to the variable in the modified Bryson-Frazier form,
/// `T = Df - Pf Psi^`, `Psi^ <- F^T (H^T S^-1 H Dp + (I - K H)^T Psi^)`.
///
/// @param states The states after the node in filter order; trailing states
///        without smoothed parameters are skipped
/// @param df The sensitivity of the filtered node to the variable
/// @param u The change of the variable
/// @param w The change of its covariance
/// @param logger The logger
/// @return Success, or `KalmanFitterError::SmoothFailed` if a state lacks
///         components or a result is not finite
template <std::size_t C>
Result<void> conditionStatesAfterNode(
    std::span<AnyMutableTrackStateProxy> states,
    const Matrix<eBoundSize, C>& df, const Vector<C>& u,
    const SquareMatrix<C>& w, const Logger& logger);

/// Condition the smoothed states before a node on its change `dx = G u`,
/// `dP = G W G^T` with the modified Bryson-Frazier adjoint:
/// `x_j += Pf_j L_j u`, `P_j += Pf_j L_j W L_j^T Pf_j`, `L <- (I - K H)^T F^T
/// L` from `L = Pp^-1 G` at the node.
///
/// @param states The states up to and including the node, in filter order
/// @param lambda `Pp^-1 G` at the node
/// @param u The change in factor form
/// @param w The covariance change in factor form
/// @param logger The logger
/// @return Success, or `KalmanFitterError::SmoothFailed` if a state lacks
///         components or a result is not finite
template <std::size_t C>
Result<void> conditionStatesBeforeNode(
    std::span<AnyMutableTrackStateProxy> states, Matrix<eBoundSize, C> lambda,
    const Vector<C>& u, const SquareMatrix<C>& w, const Logger& logger);

/// Condition the smoothed states on all measurements but the one at
/// @p rejected. Predicted and filtered states are not changed.
/// @param states The states of the track in filter order
/// @param rejected The position of the rejected state in @p states
/// @param logger The logger
/// @return The chi2 of the measurement against the smoothed state without it
Result<double> rejectSmoothedMeasurement(
    std::span<AnyMutableTrackStateProxy> states, std::size_t rejected,
    const Logger& logger);

/// Condition the smoothed states on the measurements alone: the prior of the
/// filter, the predicted state of the first state, is removed like a
/// measurement of @p parameters.
/// @param states The states of the track in filter order
/// @param parameters The parameters whose prior is removed
/// @param logger The logger
/// @return The chi2 of the prior against the smoothed state without it
Result<double> removeSmoothedPrior(std::span<AnyMutableTrackStateProxy> states,
                                   std::span<const BoundIndices> parameters,
                                   const Logger& logger);

}  // namespace Acts::detail
