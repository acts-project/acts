// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

namespace Acts {

Result<void> GainMatrixUpdater::operator()(
    const GeometryContext& gctx, MultiTrajectory::TrackStateProxy trackState,
    NavigationDirection direction, LoggerWrapper logger) const {
  (void)gctx;
  ACTS_VERBOSE("Invoked GainMatrixUpdater");

  // we should definitely have an uncalibrated measurement here
  assert(trackState.hasUncalibrated());
  // there should be a calibrated measurement
  assert(trackState.hasCalibrated());
  // we should have predicted state set
  assert(trackState.hasPredicted());
  // filtering should not have happened yet, but is allocated, therefore set
  assert(trackState.hasFiltered());

  // read-only handles. Types are eigen maps to backing storage
  const auto predicted = trackState.predicted();
  const auto predictedCovariance = trackState.predictedCovariance();

  ACTS_VERBOSE("Predicted parameters: " << predicted.transpose());
  ACTS_VERBOSE("Predicted covariance:\n" << predictedCovariance);

  // read-write handles. Types are eigen maps into backing storage.
  // This writes directly into the trajectory storage
  auto filtered = trackState.filtered();
  auto filteredCovariance = trackState.filteredCovariance();

  // default-constructed error represents success, i.e. an invalid error code
  std::error_code error;
  visit_measurement(
      trackState.calibrated(), trackState.calibratedCovariance(),
      trackState.calibratedSize(),
      [&](const auto calibrated, const auto calibratedCovariance) {
        constexpr size_t kMeasurementSize =
            decltype(calibrated)::RowsAtCompileTime;
        using ParametersVector = ActsVector<kMeasurementSize>;
        using CovarianceMatrix = ActsSymMatrix<kMeasurementSize>;

        ACTS_VERBOSE("Measurement dimension: " << kMeasurementSize);
        ACTS_VERBOSE("Calibrated measurement: " << calibrated.transpose());
        ACTS_VERBOSE("Calibrated measurement covariance:\n"
                     << calibratedCovariance);

        const auto H =
            trackState.projector()
                .template topLeftCorner<kMeasurementSize, eBoundSize>()
                .eval();

        ACTS_VERBOSE("Measurement projector H:\n" << H);

        const auto K =
            (predictedCovariance * H.transpose() *
             (H * predictedCovariance * H.transpose() + calibratedCovariance)
                 .inverse())
                .eval();

        ACTS_VERBOSE("Gain Matrix K:\n" << K);

        if (K.hasNaN()) {
          error =
              (direction == NavigationDirection::Forward)
                  ? KalmanFitterError::ForwardUpdateFailed
                  : KalmanFitterError::BackwardUpdateFailed;  // set to error
          return false;                                       // abort execution
        }

        filtered = predicted + K * (calibrated - H * predicted);
        filteredCovariance =
            (BoundSymMatrix::Identity() - K * H) * predictedCovariance;
        ACTS_VERBOSE("Filtered parameters: " << filtered.transpose());
        ACTS_VERBOSE("Filtered covariance:\n" << filteredCovariance);

        // calculate filtered residual
        //
        // FIXME: Without separate residual construction and assignment, we
        //        currently take a +0.7GB build memory consumption hit in the
        //        EventDataView unit tests. Revisit this once Measurement
        //        overhead problems (Acts issue #350) are sorted out.
        //
        ParametersVector residual;
        residual = calibrated - H * filtered;
        ACTS_VERBOSE("Residual: " << residual.transpose());

        trackState.chi2() =
            (residual.transpose() *
             ((CovarianceMatrix::Identity() - H * K) * calibratedCovariance)
                 .inverse() *
             residual)
                .value();

        ACTS_VERBOSE("Chi2: " << trackState.chi2());
        return true;  // continue execution
      });

  return error ? Result<void>::failure(error) : Result<void>::success();
}
}  // namespace Acts