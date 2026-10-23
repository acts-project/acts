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
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
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

/// @brief Condition the smoothed states of a segment on a new estimate of
/// one of its states
///
/// The segment runs from @p sharedIndex to @p tipIndex in filter order and is
/// smoothed. Each state is updated by its regression on the shared state,
/// `x_k += L_k (x_s' - x_s)`, `P_k += L_k (P_s' - P_s) L_k^T`, which is the
/// exact linear-Gaussian update for information that enters through the
/// shared state. `L_k` is the sensitivity of the smoothed state to the prior
/// mean of the segment, computed in the modified Bryson-Frazier form with a
/// forward sweep `Dp_k = F_k Df_(k-1)`, `Df_k = (I - K_k H_k) Dp_k` and a
/// backward sweep `Psi~_k = H_k^T S_k^-1 H_k Dp_k + (I - K_k H_k)^T Psi^_k`,
/// `Psi^_(k-1) = F_k^T Psi~_k`, `T_k = Df_k - Pf_k Psi^_k`, and
/// `L_k = T_k Pp_s Ps_s^-1`.
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
/// @return Success, or @ref KalmanFitterError::SmoothFailed if the shared
///         state is not on the chain, a state lacks components, or a result
///         is not finite
template <typename traj_t>
Result<void> conditionSmoothedStatesOnSharedState(
    traj_t &trajectory, std::size_t sharedIndex, std::size_t tipIndex,
    const BoundVector &newMean, const BoundMatrix &newCov,
    const Logger &logger = getDummyLogger()) {
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

  const BoundVector sharedMean = shared.smoothed();
  const BoundMatrix sharedCov = shared.smoothedCovariance();
  const BoundVector deltaMean = subtractBoundParameters(newMean, sharedMean);
  const BoundMatrix deltaCov = newCov - sharedCov;
  const Eigen::LDLT<BoundMatrix> sharedLdlt(sharedCov);
  if (sharedLdlt.info() != Eigen::Success || !sharedLdlt.isPositive()) {
    ACTS_DEBUG("Shared state covariance is not positive definite");
    return KalmanFitterError::SmoothFailed;
  }
  // Z = Pp_s Ps_s^-1, applied to the change of the shared state
  const BoundMatrix sharedPredictedCov = shared.predictedCovariance();
  const BoundMatrix z = sharedLdlt.solve(sharedPredictedCov).transpose();
  const BoundVector u = z * deltaMean;
  const BoundMatrix w = z * deltaCov * z.transpose();

  // Df = (I - K H) Dp
  auto filteredSensitivity = [](const AnyConstTrackStateProxy &ts,
                                const BoundMatrix &dp) -> BoundMatrix {
    if (!ts.typeFlags().isMeasurement()) {
      return dp;
    }
    BoundMatrix df = dp;
    visit_measurement(
        ts.calibratedSize(),
        [&]<std::size_t N>(std::integral_constant<std::size_t, N>) {
          const auto H = ts.projectorSubspaceHelper<N>().projector();
          const BoundMatrix pp = ts.predictedCovariance();
          const Eigen::Matrix<double, N, N> s =
              H * pp * H.transpose() + ts.calibratedCovariance<N>();
          const Eigen::Matrix<double, eBoundSize, N> k =
              pp * H.transpose() * s.inverse();
          df = dp - k * (H * dp);
        });
    return df;
  };
  // Psi~ = H^T S^-1 H Dp + (I - K H)^T Psi^
  auto adjointSensitivity = [](const AnyConstTrackStateProxy &ts,
                               const BoundMatrix &dp,
                               const BoundMatrix &psiHat) -> BoundMatrix {
    if (!ts.typeFlags().isMeasurement()) {
      return psiHat;
    }
    BoundMatrix psiTilde = psiHat;
    visit_measurement(
        ts.calibratedSize(),
        [&]<std::size_t N>(std::integral_constant<std::size_t, N>) {
          const auto H = ts.projectorSubspaceHelper<N>().projector();
          const BoundMatrix pp = ts.predictedCovariance();
          const Eigen::Matrix<double, N, N> s =
              H * pp * H.transpose() + ts.calibratedCovariance<N>();
          const Eigen::Matrix<double, N, N> sInv = s.inverse();
          const Eigen::Matrix<double, eBoundSize, N> k =
              pp * H.transpose() * sInv;
          psiTilde = H.transpose() * (sInv * (H * dp)) + psiHat -
                     H.transpose() * (k.transpose() * psiHat);
        });
    return psiTilde;
  };

  // forward sweep: sensitivities (Dp, Df) of the predicted and filtered
  // states to the prior mean of the segment
  std::vector<std::pair<BoundMatrix, BoundMatrix>> sensitivities;
  sensitivities.resize(indices.size());
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

  // backward sweep: sensitivity of the smoothed states and the update
  BoundMatrix psiHat = BoundMatrix::Zero();
  for (std::size_t k = 0; k < indices.size(); ++k) {
    auto ts = trajectory.getTrackState(indices[k]);
    const bool isShared = (k + 1 == indices.size());
    const auto &[dp, dfk] = sensitivities[k];

    if (!isShared) {
      const BoundMatrix filteredCov =
          ts.hasFiltered() ? BoundMatrix(ts.filteredCovariance())
                           : BoundMatrix(ts.predictedCovariance());
      const BoundMatrix t = dfk - filteredCov * psiHat;

      const BoundVector mean = normalizeBoundParameters(ts.smoothed() + t * u);
      const BoundMatrix cov = ts.smoothedCovariance() + t * w * t.transpose();
      if (!mean.allFinite() || !cov.allFinite()) {
        ACTS_DEBUG("Conditioned state " << ts.index() << " is not finite");
        return KalmanFitterError::SmoothFailed;
      }
      ts.smoothed() = mean;
      ts.smoothedCovariance() = 0.5 * (cov + cov.transpose());

      const BoundMatrix psiTilde =
          adjointSensitivity(AnyConstTrackStateProxy{ts}, dp, psiHat);
      psiHat = ts.jacobian().transpose() * psiTilde;
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

/// @brief Link a second tracklet, filtered from a state of a first tracklet,
/// behind that state
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
                                              track_proxy_t &second,
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
