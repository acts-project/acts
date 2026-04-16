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

  visit_measurement(
      ts.calibratedSize(),
      [&, this]<std::size_t N>(std::integral_constant<std::size_t, N>) {
        visitMeasurementImpl<N>(ts, bigLambdaHat, smallLambdaHat);
      });
}

}  // namespace Acts
