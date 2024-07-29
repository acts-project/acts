// This file is part of the Acts project.
//
// Copyright (C) 2021-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/detail/GsfUtils.hpp"

#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"

#include <cstddef>

namespace Acts::detail {

using TrackStateTraits =
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax, true>;

ActsScalar calculateDeterminant(
    const double* fullCalibratedCovariance,
    TrackStateTraits::Covariance predictedCovariance,
    ProjectorMapping projector, unsigned int calibratedSize) {
  return visit_measurement(calibratedSize, [&](auto N) {
    constexpr std::size_t kMeasurementSize = decltype(N)::value;

    typename Acts::TrackStateTraits<
        kMeasurementSize, true>::CalibratedCovariance calibratedCovariance{
        fullCalibratedCovariance};

    SubspaceHelper<eBoundSize> subspaceHelper(
        {projector.begin(), projector.begin() + kMeasurementSize});

    const auto H = subspaceHelper.template projector<kMeasurementSize>();

    return (H * predictedCovariance * H.transpose() + calibratedCovariance)
        .determinant();
  });
}
}  // namespace Acts::detail
