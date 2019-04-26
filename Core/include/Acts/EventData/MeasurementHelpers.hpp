// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <variant>
#include "Acts/EventData/Measurement.hpp"

namespace Acts {

class Surface;

namespace MeasurementHelpers {

/// @brief Extract surface from a type erased measurement object
/// @tparam T The FittableMeasurement type
/// @return const pointer to the extracted surface
template <typename T>
const Surface* getSurface(const T& fittable_measurement) {
  return std::visit([](const auto& meas) { return &meas.referenceSurface(); },
                    fittable_measurement);
}

template <typename T>
size_t getSize(const T& fittable_measurement) {
  return std::visit([](const auto& meas) { return meas.size(); },
                    fittable_measurement);
}
}  // namespace MeasurementHelpers
}  // namespace Acts
