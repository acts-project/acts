// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/MbfSmoother.hpp"

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/KalmanFitterError.hpp"

namespace Acts {

Result<void> MbfSmoother::calculateSmoothed(
    const GeometryContext& gctx, AnyMutableTrackStateProxy& ts,
    const BoundMatrix& bigLambdaHat, const BoundVector& smallLambdaHat) const {
  const auto filteredCovariance = ts.filteredCovariance();
  const BoundMatrix covariance = filteredCovariance - filteredCovariance *
                                                          bigLambdaHat *
                                                          filteredCovariance;
  const BoundVector parameters =
      ts.filtered() - filteredCovariance * smallLambdaHat;
  if (!parameters.allFinite() || !covariance.allFinite()) {
    return KalmanFitterError::SmoothFailed;
  }
  const auto parametrization =
      ts.hasReferenceSurface()
          ? ts.referenceSurface().boundParametrization(gctx)
          : BoundParametrization::globalSpherical();
  const auto canonical = parametrization.canonicalize(parameters, covariance);
  if (!canonical.ok()) {
    return canonical.error();
  }
  ts.smoothed() = canonical->first;
  ts.smoothedCovariance() = canonical->second;
  return Result<void>::success();
}

Result<void> MbfSmoother::visitNonMeasurement(
    const AnyConstTrackStateProxy::ConstCovarianceMap& jacobian,
    BoundMatrix& bigLambdaHat, BoundVector& smallLambdaHat) const {
  const auto F = jacobian;

  bigLambdaHat = (F.transpose() * bigLambdaHat * F).eval();
  smallLambdaHat = (F.transpose() * smallLambdaHat).eval();
  if (!bigLambdaHat.allFinite() || !smallLambdaHat.allFinite()) {
    return KalmanFitterError::SmoothFailed;
  }
  return Result<void>::success();
}

Result<void> MbfSmoother::visitMeasurement(const GeometryContext& gctx,
                                           const AnyConstTrackStateProxy& ts,
                                           BoundMatrix& bigLambdaHat,
                                           BoundVector& smallLambdaHat) const {
  assert(ts.hasCalibrated());

  return visit_measurement(
      ts.calibratedSize(),
      [&, this]<std::size_t N>(std::integral_constant<std::size_t, N>) {
        return visitMeasurementImpl<N>(gctx, ts, bigLambdaHat, smallLambdaHat);
      });
}

}  // namespace Acts
