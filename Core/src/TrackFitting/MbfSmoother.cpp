// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/MbfSmoother.hpp"

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"

#include <cstdint>

namespace Acts {

void MbfSmoother::calculateSmoothed(AnyMutableTrackStateProxy& ts,
                                    const BoundMatrix& bigLambdaHat,
                                    const BoundVector& smallLambdaHat) const {
  auto filteredCovariance = ts.filteredCovariance();
  auto smoothed = ts.smoothed();
  ts.smoothedCovariance() = filteredCovariance - filteredCovariance *
                                                     bigLambdaHat *
                                                     filteredCovariance;
  smoothed = ts.filtered() - filteredCovariance * smallLambdaHat;
  // Normalize phi and theta
  smoothed = normalizeBoundParameters(smoothed);
}

void MbfSmoother::visitNonMeasurement(
    const AnyConstTrackStateProxy::ConstCovarianceMap& jacobian,
    BoundMatrix& bigLambdaHat, BoundVector& smallLambdaHat) const {
  const auto F = jacobian;

  bigLambdaHat = F.transpose() * bigLambdaHat * F;
  smallLambdaHat = F.transpose() * smallLambdaHat;
}

void MbfSmoother::visitMeasurement(const AnyConstTrackStateProxy& ts,
                                   BoundMatrix& bigLambdaHat,
                                   BoundVector& smallLambdaHat) const {
  assert(ts.hasCalibrated());

  const auto F = ts.jacobian();

  visit_measurement(
      ts.calibratedSize(),
      [&]<std::size_t N>(std::integral_constant<std::size_t, N> /*unused*/) {
        constexpr std::size_t kMeasurementSize = N;

        const auto subspaceHelper =
            ts.projectorSubspaceHelper<kMeasurementSize>();

        using ProjectorMatrix =
            Eigen::Matrix<double, kMeasurementSize, eBoundSize>;
        using CovarianceMatrix =
            Eigen::Matrix<double, kMeasurementSize, kMeasurementSize>;
        using KalmanGainMatrix =
            Eigen::Matrix<double, eBoundSize, kMeasurementSize>;

        typename TrackStateTraits<kMeasurementSize, true>::Calibrated
            calibrated{ts.calibrated<kMeasurementSize>()};
        typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
            calibratedCovariance{ts.calibratedCovariance<kMeasurementSize>()};

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
        const Eigen::Matrix<double, kMeasurementSize, 1> y =
            (calibrated - H * ts.predicted());

        const Acts::BoundMatrix bigLambdaTilde =
            (H.transpose() * SInv * H + CHat.transpose() * bigLambdaHat * CHat);
        const Eigen::Matrix<double, eBoundSize, 1> smallLambdaTilde =
            (-H.transpose() * SInv * y + CHat.transpose() * smallLambdaHat);

        bigLambdaHat = F.transpose() * bigLambdaTilde * F;
        smallLambdaHat = F.transpose() * smallLambdaTilde;
      });
}

}  // namespace Acts
