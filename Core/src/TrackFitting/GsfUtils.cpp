// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/GsfUtils.hpp"

#include "Acts/EventData/MeasurementHelpers.hpp"

namespace Acts {
namespace detail {

using TrackStateTraits =
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax, true>;

ActsScalar calculateDeterminant(
    TrackStateTraits::Measurement fullCalibrated,
    TrackStateTraits::MeasurementCovariance fullCalibratedCovariance,
    TrackStateTraits::Covariance predictedCovariance,
    TrackStateTraits::Projector projector, unsigned int calibratedSize) {
  return visit_measurement(
      fullCalibrated, fullCalibratedCovariance, calibratedSize,
      [&](const auto calibrated, const auto calibratedCovariance) {
        constexpr size_t kMeasurementSize =
            decltype(calibrated)::RowsAtCompileTime;
        const auto H =
            projector.template topLeftCorner<kMeasurementSize, eBoundSize>()
                .eval();

        return (H * predictedCovariance * H.transpose() + calibratedCovariance)
            .determinant();
      });
}
}  // namespace detail
}  // namespace Acts
