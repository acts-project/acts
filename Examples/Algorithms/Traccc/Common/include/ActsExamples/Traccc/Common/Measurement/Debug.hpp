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
#include "ActsExamples/Traccc/Common/Conversion/MeasurementMatch.hpp"
#include "ActsExamples/Traccc/Common/Util/IndexMap.hpp"

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

// This file is for debugging and for getting the matching between two
// collections of measurements as a string.

namespace ActsExamples::Traccc::Common::Measurement {

/// @brief Creates a string with the data of the measurements and their relation according to the index map.
/// @param measurements1 the measurements (1).
/// @param measurements2 the measurements (2).
/// @param indexMap the index map: measurements1 indices -> measurement2 indices.
/// The index map describes which elements are related in the two measurement
/// collections.
/// @return a string formatted as a table.
std::string pairingStatistics(
    const std::vector<ActsExamples::BoundVariantMeasurement>& measurements1,
    const std::vector<ActsExamples::BoundVariantMeasurement>& measurements2,
    const std::map<std::size_t, std::size_t>& indexMap);

template <typename allocator_t, typename detector_t>
std::string pairingStatistics(
    const std::vector<traccc::measurement, allocator_t>& tracccMeasurements,
    const std::vector<ActsExamples::BoundVariantMeasurement>& measurements,
    const detector_t& detector)
    {

    std::vector<ActsExamples::BoundVariantMeasurement> convertedMeasurements;
    Conversion::convertMeasurements(detector, tracccMeasurements, convertedMeasurements);
    
    auto indexMap = Util::matchMap(convertedMeasurements, measurements, Conversion::MeasurementGeoIDHash{}, Conversion::MeasurementAproxEquals<>{});

    return pairingStatistics(convertedMeasurements,
                             measurements, indexMap);
  }

}  // namespace ActsExamples::Traccc::Common::Measurement
