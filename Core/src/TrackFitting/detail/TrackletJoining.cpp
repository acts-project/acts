// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/TrackletJoining.hpp"

#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"

#include <Eigen/Cholesky>

namespace Acts::detail {

Result<std::pair<BoundVector, BoundMatrix>> smoothPriorNode(
    const AnyConstTrackStateProxy& first, const BoundVector& priorMean,
    const BoundMatrix& priorCov, const Logger& logger) {
  const BoundMatrix predictedCov = first.predictedCovariance();
  const Eigen::LDLT<BoundMatrix> ldlt(predictedCov);
  if (ldlt.info() != Eigen::Success || !ldlt.isPositive()) {
    ACTS_DEBUG("Predicted covariance of state " << first.index()
                                                << " is not positive definite");
    return KalmanFitterError::SmoothFailed;
  }

  // G = P F^T Pp^-1
  const BoundMatrix jacobian = first.jacobian();
  const BoundMatrix g = ldlt.solve(jacobian * priorCov).transpose();

  const BoundVector mean = normalizeBoundParameters(
      priorMean +
      g * subtractBoundParameters(first.smoothed(), first.predicted()));
  const BoundMatrix cov =
      priorCov + g * (BoundMatrix(first.smoothedCovariance()) - predictedCov) *
                     g.transpose();
  if (!mean.allFinite() || !cov.allFinite()) {
    ACTS_DEBUG("Prior node estimate is not finite");
    return KalmanFitterError::SmoothFailed;
  }
  return std::pair{mean, BoundMatrix(0.5 * (cov + cov.transpose()))};
}

Result<SharedStateCondition> makeSharedStateCondition(
    const AnyConstTrackStateProxy& shared, const BoundVector& newMean,
    const BoundMatrix& newCov, const Logger& logger) {
  const BoundMatrix smoothedCov = shared.smoothedCovariance();
  const Eigen::LDLT<BoundMatrix> ldlt(smoothedCov);
  if (ldlt.info() != Eigen::Success || !ldlt.isPositive()) {
    ACTS_DEBUG("Shared state covariance is not positive definite");
    return KalmanFitterError::SmoothFailed;
  }

  // Z = Pp Ps^-1, applied to the change of the shared state
  const BoundMatrix predictedCov = shared.predictedCovariance();
  const BoundMatrix z = ldlt.solve(predictedCov).transpose();
  const BoundVector deltaMean =
      subtractBoundParameters(newMean, shared.smoothed());
  const BoundMatrix deltaCov = newCov - smoothedCov;

  SharedStateCondition condition;
  condition.u = z * deltaMean;
  condition.w = z * deltaCov * z.transpose();
  return condition;
}

BoundMatrix filteredSensitivity(const AnyConstTrackStateProxy& ts,
                                const BoundMatrix& dp) {
  if (!ts.typeFlags().isMeasurement()) {
    return dp;
  }
  BoundMatrix df = dp;
  visit_measurement(
      ts.calibratedSize(),
      [&]<std::size_t N>(std::integral_constant<std::size_t, N>) {
        const auto H = ts.projectorSubspaceHelper<N>().projector();
        const BoundMatrix pp = ts.predictedCovariance();
        const SquareMatrix<N> s =
            H * pp * H.transpose() + ts.calibratedCovariance<N>();
        const Matrix<eBoundSize, N> k = pp * H.transpose() * s.inverse();
        // Df = (I - K H) Dp
        df = dp - k * (H * dp);
      });
  return df;
}

Result<void> conditionTrackState(AnyMutableTrackStateProxy ts,
                                 const BoundMatrix& dp, const BoundMatrix& df,
                                 BoundMatrix& psiHat,
                                 const SharedStateCondition& condition,
                                 const Logger& logger) {
  const BoundMatrix filteredCov = ts.hasFiltered()
                                      ? BoundMatrix(ts.filteredCovariance())
                                      : BoundMatrix(ts.predictedCovariance());
  // T = Df - Pf Psi^
  const BoundMatrix t = df - filteredCov * psiHat;

  const BoundVector mean =
      normalizeBoundParameters(ts.smoothed() + t * condition.u);
  const BoundMatrix cov =
      ts.smoothedCovariance() + t * condition.w * t.transpose();
  if (!mean.allFinite() || !cov.allFinite()) {
    ACTS_DEBUG("Conditioned state " << ts.index() << " is not finite");
    return KalmanFitterError::SmoothFailed;
  }
  ts.smoothed() = mean;
  ts.smoothedCovariance() = 0.5 * (cov + cov.transpose());

  // Psi~ = H^T S^-1 H Dp + (I - K H)^T Psi^, then transport to the previous
  BoundMatrix psiTilde = psiHat;
  if (ts.typeFlags().isMeasurement()) {
    visit_measurement(
        ts.calibratedSize(),
        [&]<std::size_t N>(std::integral_constant<std::size_t, N>) {
          const auto H = ts.projectorSubspaceHelper<N>().projector();
          const BoundMatrix pp = ts.predictedCovariance();
          const SquareMatrix<N> s =
              H * pp * H.transpose() + ts.calibratedCovariance<N>();
          const SquareMatrix<N> sInv = s.inverse();
          const Matrix<eBoundSize, N> k = pp * H.transpose() * sInv;
          psiTilde = H.transpose() * (sInv * (H * dp)) + psiHat -
                     H.transpose() * (k.transpose() * psiHat);
        });
  }
  psiHat = ts.jacobian().transpose() * psiTilde;
  return Result<void>::success();
}

}  // namespace Acts::detail
