// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

namespace Acts {

std::tuple<double, std::error_code> GainMatrixUpdater::visitMeasurement(
    InternalTrackState trackState, NavigationDirection direction,
    LoggerWrapper logger) const {
  // default-constructed error represents success, i.e. an invalid error code
  std::error_code error;
  double chi2 = 0;
  visit_measurement(
      trackState.calibrated, trackState.calibratedCovariance,
      trackState.calibratedSize,
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
            trackState.projector
                .template topLeftCorner<kMeasurementSize, eBoundSize>()
                .eval();

        ACTS_VERBOSE("Measurement projector H:\n" << H);

        const auto K = (trackState.predictedCovariance * H.transpose() *
                        (H * trackState.predictedCovariance * H.transpose() +
                         calibratedCovariance)
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

        trackState.filtered =
            trackState.predicted + K * (calibrated - H * trackState.predicted);
        trackState.filteredCovariance = (BoundSymMatrix::Identity() - K * H) *
                                        trackState.predictedCovariance;
        ACTS_VERBOSE(
            "Filtered parameters: " << trackState.filtered.transpose());
        ACTS_VERBOSE("Filtered covariance:\n" << trackState.filteredCovariance);

        // calculate filtered residual
        //
        // FIXME: Without separate residual construction and assignment, we
        //        currently take a +0.7GB build memory consumption hit in the
        //        EventDataView unit tests. Revisit this once Measurement
        //        overhead problems (Acts issue #350) are sorted out.
        //
        ParametersVector residual;
        residual = calibrated - H * trackState.filtered;
        ACTS_VERBOSE("Residual: " << residual.transpose());

        chi2 = (residual.transpose() *
                ((CovarianceMatrix::Identity() - H * K) * calibratedCovariance)
                    .inverse() *
                residual)
                   .value();

        ACTS_VERBOSE("Chi2: " << chi2);
        return true;  // continue execution
      });

  return {chi2, error};
}

}  // namespace Acts
