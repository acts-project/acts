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
#include "ActsExamples/Traccc/Common/Conversion/MeasurementConversion.hpp"
#include "ActsExamples/Traccc/Common/Util/IndexMap.hpp"
#include "ActsExamples/Traccc/Common/Util/MapUtil.hpp"

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
template <double max_dist = .001>
struct MeasurementAproxEquals {

    bool operator()(const ActsExamples::BoundVariantMeasurement& measurement1,
    const ActsExamples::BoundVariantMeasurement& measurement2) const {
        auto gidEq = Conversion::getGeometryID(measurement1) == Conversion::getGeometryID(measurement2);

        auto sqNorm =
            (Conversion::getLocal(measurement1) - Conversion::getLocal(measurement2))
                .squaredNorm();
        auto locEq = sqNorm <= max_dist * max_dist;

        return gidEq && locEq;
    }
};


/// @brief Generates a hash for the measurement.
/// This hash is used for the locality sensitive hashing to calculate the match
/// map. Thus, this hash is not sensitive to small variations in position that
/// could could from numerical errors.
struct MeasurementGeoIDHash {
    bool operator()(const ActsExamples::BoundVariantMeasurement& measurement) const {
        return static_cast<std::size_t>(Conversion::getGeometryID(measurement).value());
    }
};

/// @brief Creates a map from traccc measurements to acts measurements.
/// The traccc elements will map to equivalent acts measurements.
/// The resulting map is assumed to be bijective, thus if any element is
/// unable to find a match an error is thrown.
/// The map is created by matching the elements in two collection based on their geometry id and approximate position.
/// The two collections of measurements are assumed to contain the same
/// measurements. However, equivalent measurements may be at different indices
/// in the collections. This function determines which indexes correspond to
// matching measurements.
/// @param tracccMeasurements the traccc measurements.
/// @param measurements the acts measurements.
template <typename allocator_t, typename detector_t>
auto matchMeasurements(
    const std::vector<traccc::measurement, allocator_t>& tracccMeasurements,
    const std::vector<ActsExamples::BoundVariantMeasurement>& measurements,
    const detector_t& detector)
    {

  std::vector<ActsExamples::BoundVariantMeasurement> convertedMeasurements;
  Conversion::convertMeasurements(detector, tracccMeasurements, convertedMeasurements);

  auto indexMap = Util::matchMap(convertedMeasurements, measurements, Conversion::MeasurementGeoIDHash{}, Conversion::MeasurementAproxEquals<>{});

  return Util::createFromIndexMap<TracccMeasurementHash, std::equal_to<traccc::measurement>>(tracccMeasurements, measurements, indexMap);
}

}  // namespace ActsExamples::Traccc::Common::Conversion
