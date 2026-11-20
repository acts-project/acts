// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/TrackletJoining.hpp"

#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"

#include <Eigen/Cholesky>

namespace {

using namespace Acts;

/// Inverse residual covariance and Kalman gain of a measurement state
template <std::size_t N>
struct MeasurementGain {
  Matrix<N, eBoundSize> h;
  SquareMatrix<N> sInv;
  Matrix<eBoundSize, N> k;

  template <typename track_state_proxy_t>
  explicit MeasurementGain(const track_state_proxy_t& ts)
      : h(ts.template projectorSubspaceHelper<N>().projector()) {
    const BoundMatrix pp = ts.predictedCovariance();
    const SquareMatrix<N> s =
        h * pp * h.transpose() + ts.template calibratedCovariance<N>();
    sInv = s.inverse();
    k = pp * h.transpose() * sInv;
  }
};

}  // namespace

namespace Acts::detail {

void setSmoothedWithoutSharing(AnyMutableTrackStateProxy ts, BoundVector mean,
                               BoundMatrix cov) {
  if (ts.hasSmoothed() &&
      ((ts.hasFiltered() && ts.smoothed().data() == ts.filtered().data()) ||
       (ts.hasPredicted() && ts.smoothed().data() == ts.predicted().data()))) {
    ts.unset(TrackStatePropMask::Smoothed);
  }
  ts.addComponents(TrackStatePropMask::Smoothed);
  ts.smoothed() = mean;
  ts.smoothedCovariance() = cov;
}

Result<std::pair<BoundVector, BoundMatrix>> smoothPriorNode(
    const GeometryContext& gctx, const AnyConstTrackStateProxy& first,
    const BoundParametrization& priorParametrization,
    const BoundVector& priorMean, const BoundMatrix& priorCov,
    const Logger& logger) {
  const BoundMatrix predictedCov = first.predictedCovariance();
  if (!predictedCov.allFinite()) {
    return KalmanFitterError::SmoothFailed;
  }
  const Eigen::LDLT<BoundMatrix> ldlt(predictedCov);
  if (ldlt.info() != Eigen::Success || !(ldlt.vectorD().array() > 0).all()) {
    ACTS_DEBUG("Predicted covariance of state " << first.index()
                                                << " is not positive definite");
    return KalmanFitterError::SmoothFailed;
  }

  // G = P F^T Pp^-1
  const BoundMatrix jacobian = first.jacobian();
  const BoundMatrix g = ldlt.solve(jacobian * priorCov).transpose();

  // The residual lives in the basis of the first state; G maps it to the prior
  const BoundVector mean =
      priorMean + g * trackletStateParametrization(gctx, first)
                          .difference(first.smoothed(), first.predicted());
  const BoundMatrix cov =
      priorCov + g * (BoundMatrix(first.smoothedCovariance()) - predictedCov) *
                     g.transpose();
  if (!mean.allFinite() || !cov.allFinite()) {
    ACTS_DEBUG("Prior node estimate is not finite");
    return KalmanFitterError::SmoothFailed;
  }
  return priorParametrization.canonicalize(
      mean, BoundMatrix(0.5 * (cov + cov.transpose())));
}

Result<SharedStateCondition> makeSharedStateCondition(
    const GeometryContext& gctx, const AnyConstTrackStateProxy& shared,
    const BoundVector& newMean, const BoundMatrix& newCov,
    const Logger& logger) {
  const BoundMatrix smoothedCov = shared.smoothedCovariance();
  if (!smoothedCov.allFinite()) {
    return KalmanFitterError::SmoothFailed;
  }
  const Eigen::LDLT<BoundMatrix> ldlt(smoothedCov);
  if (ldlt.info() != Eigen::Success || !(ldlt.vectorD().array() > 0).all()) {
    ACTS_DEBUG("Shared state covariance is not positive definite");
    return KalmanFitterError::SmoothFailed;
  }

  // Z = Pp Ps^-1, applied to the change of the shared state
  const BoundMatrix predictedCov = shared.predictedCovariance();
  const BoundMatrix z = ldlt.solve(predictedCov).transpose();
  const BoundVector deltaMean = trackletStateParametrization(gctx, shared)
                                    .difference(newMean, shared.smoothed());
  const BoundMatrix deltaCov = newCov - smoothedCov;

  SharedStateCondition condition;
  condition.u = z * deltaMean;
  condition.w = z * deltaCov * z.transpose();
  if (!condition.u.allFinite() || !condition.w.allFinite()) {
    ACTS_DEBUG("Shared state conditioning terms are not finite");
    return KalmanFitterError::SmoothFailed;
  }
  return condition;
}

BoundMatrix filteredSensitivity(const AnyConstTrackStateProxy& ts,
                                const BoundMatrix& dp) {
  if (!ts.typeFlags().isMeasurement()) {
    return dp;
  }
  BoundMatrix df = dp;
  visit_measurement(ts.calibratedSize(),
                    [&]<std::size_t N>(std::integral_constant<std::size_t, N>) {
                      const MeasurementGain<N> gain(ts);
                      df = dp - gain.k * (gain.h * dp);
                    });
  return df;
}

Result<void> conditionTrackState(const GeometryContext& gctx,
                                 AnyMutableTrackStateProxy ts,
                                 const BoundMatrix& dp, const BoundMatrix& df,
                                 BoundMatrix& psiHat,
                                 const SharedStateCondition& condition,
                                 const Logger& logger) {
  const BoundMatrix filteredCov = ts.hasFiltered()
                                      ? BoundMatrix(ts.filteredCovariance())
                                      : BoundMatrix(ts.predictedCovariance());
  // T = Df - Pf Psi^
  const BoundMatrix t = df - filteredCov * psiHat;

  const BoundVector mean = ts.smoothed() + t * condition.u;
  const BoundMatrix cov =
      ts.smoothedCovariance() + t * condition.w * t.transpose();
  if (!mean.allFinite() || !cov.allFinite()) {
    ACTS_DEBUG("Conditioned state " << ts.index() << " is not finite");
    return KalmanFitterError::SmoothFailed;
  }
  const auto canonical = trackletStateParametrization(gctx, ts).canonicalize(
      mean, BoundMatrix(0.5 * (cov + cov.transpose())));
  if (!canonical.ok()) {
    return canonical.error();
  }
  setSmoothedWithoutSharing(ts, canonical->first, canonical->second);

  // Psi~ = H^T S^-1 H Dp + (I - K H)^T Psi^, then transport to the previous
  BoundMatrix psiTilde = psiHat;
  if (ts.typeFlags().isMeasurement()) {
    visit_measurement(
        ts.calibratedSize(),
        [&]<std::size_t N>(std::integral_constant<std::size_t, N>) {
          const MeasurementGain<N> gain(ts);
          psiTilde = gain.h.transpose() * (gain.sInv * (gain.h * dp)) + psiHat -
                     gain.h.transpose() * (gain.k.transpose() * psiHat);
        });
  }
  psiHat = ts.jacobian().transpose() * psiTilde;
  return Result<void>::success();
}

}  // namespace Acts::detail
