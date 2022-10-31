// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/MeasurementSelector.hpp"

namespace Acts {

double MeasurementSelector::calculateChi2(
    double* fullCalibrated, double* fullCalibratedCovariance,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Parameters predicted,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance predictedCovariance,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Projector projector,
    unsigned int calibratedSize) const {
  return visit_measurement(calibratedSize, [&](auto N) -> double {
    constexpr size_t kMeasurementSize = decltype(N)::value;

    typename TrackStateTraits<kMeasurementSize, true>::Measurement calibrated{
        fullCalibrated};

    typename TrackStateTraits<kMeasurementSize, true>::MeasurementCovariance
        calibratedCovariance{fullCalibratedCovariance};

    using ParametersVector = ActsVector<kMeasurementSize>;

    // Take the projector (measurement mapping function)
    const auto H =
        projector.template topLeftCorner<kMeasurementSize, eBoundSize>().eval();

    // Get the residuals
    ParametersVector res;
    res = calibrated - H * predicted;

    // Get the chi2
    return (res.transpose() *
            ((calibratedCovariance + H * predictedCovariance * H.transpose()))
                .inverse() *
            res)
        .eval()(0, 0);
  });
}

}  // namespace Acts
