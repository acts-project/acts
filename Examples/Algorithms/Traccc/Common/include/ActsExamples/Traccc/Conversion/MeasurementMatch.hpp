// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s)
#include "ActsExamples/EventData/Measurement.hpp"

// Acts Examples include(s)
#include "ActsExamples/Traccc/Conversion/MeasurementConversion.hpp"
#include "ActsExamples/Traccc/Util/IndexMap.hpp"
#include "ActsExamples/Traccc/Util/MapUtil.hpp"

// Traccc include(s)
#include "traccc/edm/measurement.hpp"

// System include(s)
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <variant>
#include <vector>

namespace ActsExamples::Traccc::Common::Conversion {

/// @brief Checks if two measurements are equal based on their geometry ID and position alone.
/// @note max_dist is the tolerance of distance in local position.
/// The distance between the local positions of the measurements must be
/// less or equal to this value to be considered equal.
/// @returns true or false depending on whether they are considered equal.
struct MeasurementAproxEquals {
  bool operator()(const ActsExamples::MeasurementContainer::ConstVariableProxy&
                      measurement1,
                  const ActsExamples::MeasurementContainer::ConstVariableProxy&
                      measurement2) const {
    auto gidEq = Conversion::getGeometryID(measurement1) ==
                 Conversion::getGeometryID(measurement2);

    auto sqNorm = (Conversion::getLocal(measurement1) -
                   Conversion::getLocal(measurement2))
                      .squaredNorm();
    auto locEq = sqNorm <= .001 * .001;

    return gidEq && locEq;
  }
};

/// @brief Generates a hash for the measurement.
/// This hash is used for the locality sensitive hashing to calculate the match
/// map. Thus, this hash is not sensitive to small variations in position that
/// could could from numerical errors.
struct MeasurementGeoIDHash {
  std::size_t operator()(
      const ActsExamples::MeasurementContainer::ConstVariableProxy& measurement)
      const {
    return Conversion::getGeometryID(measurement).value();
  }
};
}  // namespace ActsExamples::Traccc::Common::Conversion
