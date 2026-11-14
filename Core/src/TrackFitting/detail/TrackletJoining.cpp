// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/TrackletJoining.hpp"

#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"

#include <array>

#include <Eigen/Cholesky>

namespace Acts::detail {

Result<std::pair<BoundVector, BoundMatrix>> smoothPriorNode(
    const AnyConstTrackStateProxy& first, const BoundVector& priorMean,
    const BoundMatrix& priorCov, const Logger& logger) {
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
  return makeNodeCondition(shared.predictedCovariance(), shared.smoothed(),
                           shared.smoothedCovariance(), newMean, newCov,
                           logger);
}

Result<SharedStateCondition> makeNodeCondition(const BoundMatrix& predictedCov,
                                               const BoundVector& smoothedMean,
                                               const BoundMatrix& smoothedCov,
                                               const BoundVector& newMean,
                                               const BoundMatrix& newCov,
                                               const Logger& logger) {
  if (!smoothedCov.allFinite()) {
    return KalmanFitterError::SmoothFailed;
  }
  const Eigen::LDLT<BoundMatrix> ldlt(smoothedCov);
  if (ldlt.info() != Eigen::Success || !(ldlt.vectorD().array() > 0).all()) {
    ACTS_DEBUG("Shared state covariance is not positive definite");
    return KalmanFitterError::SmoothFailed;
  }

  // Z = Pp Ps^-1, applied to the change of the shared state
  const BoundMatrix z = ldlt.solve(predictedCov).transpose();
  const BoundVector deltaMean = subtractBoundParameters(newMean, smoothedMean);
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

Result<std::pair<BoundVector, BoundMatrix>> subtractSecondPrior(
    const BoundVector& sharedMean, const BoundMatrix& sharedCov,
    const BoundVector& priorMean, const BoundMatrix& priorCov,
    const BoundVector& nodeMean, const BoundMatrix& nodeCov,
    const Logger& logger) {
  const std::array<const BoundMatrix*, 3> covs = {&sharedCov, &priorCov,
                                                  &nodeCov};
  std::array<BoundMatrix, 3> information;
  for (std::size_t i = 0; i < covs.size(); ++i) {
    if (!covs[i]->allFinite()) {
      return KalmanFitterError::SmoothFailed;
    }
    const Eigen::LDLT<BoundMatrix> ldlt(*covs[i]);
    if (ldlt.info() != Eigen::Success || !(ldlt.vectorD().array() > 0).all()) {
      ACTS_DEBUG("Covariance at the shared state is not positive definite");
      return KalmanFitterError::SmoothFailed;
    }
    information[i] = ldlt.solve(BoundMatrix::Identity());
  }
  const auto& [yShared, yPrior, yNode] = information;

  const BoundMatrix y = yShared - yPrior + yNode;
  const Eigen::LDLT<BoundMatrix> ldlt(y);
  if (ldlt.info() != Eigen::Success || !(ldlt.vectorD().array() > 0).all()) {
    ACTS_DEBUG("Combined information is not positive definite");
    return KalmanFitterError::SmoothFailed;
  }
  const BoundMatrix cov = ldlt.solve(BoundMatrix::Identity());
  // relative to the shared mean, so that the angles need no unwrapping
  const BoundVector mean = normalizeBoundParameters(
      sharedMean +
      cov * (yNode * subtractBoundParameters(nodeMean, sharedMean) -
             yPrior * subtractBoundParameters(priorMean, sharedMean)));
  if (!mean.allFinite() || !cov.allFinite()) {
    ACTS_DEBUG("Combined estimate is not finite");
    return KalmanFitterError::SmoothFailed;
  }
  return std::pair{mean, BoundMatrix(0.5 * (cov + cov.transpose()))};
}

Result<void> conditionSmoothedStatesOnPriorNode(
    std::span<AnyMutableTrackStateProxy> states, const BoundMatrix& priorCov,
    const BoundVector& nodeMean, const BoundMatrix& nodeCov,
    const BoundVector& newMean, const BoundMatrix& newCov,
    const Logger& logger) {
  auto condition =
      makeNodeCondition(priorCov, nodeMean, nodeCov, newMean, newCov, logger);
  if (!condition.ok()) {
    return condition.error();
  }
  return conditionStatesAfterNode<eBoundSize>(
      states, BoundMatrix::Identity(), condition->u, condition->w, logger);
}

}  // namespace Acts::detail
