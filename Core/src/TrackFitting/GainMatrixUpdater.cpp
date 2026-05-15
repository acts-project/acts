// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"

#include <cstddef>
#include <type_traits>

namespace Acts {

Result<void> GainMatrixUpdater::visitMeasurement(
    AnyMutableTrackStateProxy trackState, const Logger& logger) const {
  // default-constructed error represents success, i.e. an invalid error code

  return visit_measurement(
      trackState.calibratedSize(),
      [&, this]<std::size_t N>(
          std::integral_constant<std::size_t, N>) -> Result<void> {
        return visitMeasurementImpl<N>(trackState, logger);
      });
}

}  // namespace Acts
