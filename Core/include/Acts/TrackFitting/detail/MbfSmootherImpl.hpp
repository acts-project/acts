// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/TrackFitting/MbfSmoother.hpp"

namespace Acts {

template <std::size_t N>
void MbfSmoother::visitMeasurementImpl(const AnyConstTrackStateProxy& ts,
                                       BoundMatrix& bigLambdaHat,
                                       BoundVector& smallLambdaHat) const {
  const auto F = ts.jacobian();

  const auto subspaceHelper = ts.projectorSubspaceHelper<N>();

  using ProjectorMatrix = Eigen::Matrix<double, N, eBoundSize>;
  using CovarianceMatrix = Eigen::Matrix<double, N, N>;
  using KalmanGainMatrix = Eigen::Matrix<double, eBoundSize, N>;

  typename TrackStateTraits<N, true>::Calibrated calibrated{ts.calibrated<N>()};
  typename TrackStateTraits<N, true>::CalibratedCovariance calibratedCovariance{
      ts.calibratedCovariance<N>()};

  // Projector matrix
  const ProjectorMatrix H = subspaceHelper.projector();

  // Predicted parameter covariance
  const auto predictedCovariance = ts.predictedCovariance();

  // Residual covariance
  const CovarianceMatrix S =
      (H * predictedCovariance * H.transpose() + calibratedCovariance);
  // TODO Sinv could be cached by the filter step
  const CovarianceMatrix SInv = S.inverse();

  // Kalman gain
  // TODO K could be cached by the filter step
  const KalmanGainMatrix K = (predictedCovariance * H.transpose() * SInv);

  const Acts::BoundMatrix CHat = (Acts::BoundMatrix::Identity() - K * H);
  const Eigen::Matrix<double, N, 1> y = (calibrated - H * ts.predicted());

  const Acts::BoundMatrix bigLambdaTilde =
      (H.transpose() * SInv * H + CHat.transpose() * bigLambdaHat * CHat);
  const Eigen::Matrix<double, eBoundSize, 1> smallLambdaTilde =
      (-H.transpose() * SInv * y + CHat.transpose() * smallLambdaHat);

  bigLambdaHat = F.transpose() * bigLambdaTilde * F;
  smallLambdaHat = F.transpose() * smallLambdaTilde;
}

}  // namespace Acts
