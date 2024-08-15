// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/MbfSmoother.hpp"

namespace Acts {

void MbfSmoother::calculateSmoothed(InternalTrackState& ts,
                                    const BoundMatrix& big_lambda_hat,
                                    const BoundVector& small_lambda_hat) const {
  ts.smoothedCovariance = ts.filteredCovariance - ts.filteredCovariance *
                                                      big_lambda_hat *
                                                      ts.filteredCovariance;
  ts.smoothed = ts.filtered - ts.filteredCovariance * small_lambda_hat;
}

void MbfSmoother::visitNonMeasurement(const InternalTrackState& ts,
                                      BoundMatrix& big_lambda_hat,
                                      BoundVector& small_lambda_hat) const {
  const auto F = ts.jacobian;

  big_lambda_hat = F.transpose() * big_lambda_hat * F;
  small_lambda_hat = F.transpose() * small_lambda_hat;
}

void MbfSmoother::visitMeasurement(const InternalTrackState& ts,
                                   BoundMatrix& big_lambda_hat,
                                   BoundVector& small_lambda_hat) const {
  assert(ts.measurement.has_value());

  const auto& measurement = ts.measurement.value();

  visit_measurement(measurement.calibratedSize, [&](auto N) -> void {
    constexpr std::size_t kMeasurementSize = decltype(N)::value;

    typename TrackStateTraits<kMeasurementSize, true>::Calibrated calibrated{
        measurement.calibrated};
    typename TrackStateTraits<kMeasurementSize, true>::CalibratedCovariance
        calibratedCovariance{measurement.calibratedCovariance};

    // Measurement matrix
    const auto H = measurement.projector
                       .template topLeftCorner<kMeasurementSize, eBoundSize>()
                       .eval();

    // Residual covariance
    const auto S =
        (H * ts.predictedCovariance * H.transpose() + calibratedCovariance)
            .eval();
    // TODO Sinv could be cached by the filter step
    const auto S_inv = S.inverse().eval();

    // Kalman gain
    // TODO K could be cached by the filter step
    const auto K = (ts.predictedCovariance * H.transpose() * S_inv).eval();

    const auto C_hat = (BoundMatrix::Identity() - K * H).eval();
    const auto y = (calibrated - H * ts.predicted).eval();

    const auto big_lambda_tilde =
        (H.transpose() * S_inv * H + C_hat.transpose() * big_lambda_hat * C_hat)
            .eval();
    const auto small_lambda_tilde =
        (-H.transpose() * S_inv * y + C_hat.transpose() * small_lambda_hat)
            .eval();

    const auto F = ts.jacobian;

    big_lambda_hat = F.transpose() * big_lambda_tilde * F;
    small_lambda_hat = F.transpose() * small_lambda_tilde;
  });
}

}  // namespace Acts
