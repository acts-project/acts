// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Traccc/Conversion/MeasurementConversion.hpp"

#include <cstdint>
#include <cstdlib>
#include <map>
#include <utility>
#include <variant>

#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/edm/spacepoint.hpp"

namespace ActsExamples::Traccc::Common::Conversion {

/// @brief Converts a collection of traccc space points and appends the result to the given outputContainer.
/// @param spacepoints the traccc spacepoints.
/// @param meausrementConv the measurement ConversionData (traccc measurement  -> acts bound variant measurement).
/// @param outputContainer the container to put the converted space points into (an empty container is expected).
/// @returns The spacepoint ConversionData (traccc space point -> acts space point)
std::pair<SimSpacePointContainer, std::map<std::size_t, std::size_t>>
convertTracccToActsSpacePoints(
    std::vector<traccc::spacepoint,
                std::pmr::polymorphic_allocator<traccc::spacepoint>>&
        spacePoints,
    const std::map<std::size_t, std::size_t>& tracccToActsMeasurementIndexMap,
    const MeasurementContainer& actsMeasurements);

/// @brief Converts a collection of space points to traccc space points and appends the result to the given outputContainer.
/// @param spacepoints the traccc space points.
/// @param meausrementConv the measurement ConversionData (acts bound variant measurement -> traccc measurement).
/// @param outputContainer the container to put the converted space points into (an empty container is expected).
/// @returns The spacepoint ConversionData (acts space point -> traccc space point)
std::pair<std::vector<traccc::spacepoint,
                      std::pmr::polymorphic_allocator<traccc::spacepoint>>,
          std::map<std::size_t, std::size_t>>
convertActsToTracccSpacePoints(
    const SimSpacePointContainer& actsSpacePoints,
    const std::map<std::size_t, std::size_t>& actsToTracccMeasurementIndexMap,
    const std::vector<traccc::measurement,
                      std::pmr::polymorphic_allocator<traccc::measurement>>&
        tracccMeasurements);

}  // namespace ActsExamples::Traccc::Common::Conversion
