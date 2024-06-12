// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Plugin include(s)
#include "Acts/Plugins/Traccc/MeasurementConversion.hpp"

// Acts include(s)
#include "Acts/Geometry/GeometryIdentifier.hpp"

// Acts examples include(s)
#include "ActsExamples/EventData/IndexSourceLink.hpp"

// Traccc include(s)
#include "traccc/edm/measurement.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <vector>

namespace ActsExamples::Traccc::Common::Conversion {

/// @brief Converts traccc measurements to acts measurements.
/// @param detector The detray detector,
/// @param measurements The traccc measurements,
/// @return A vector of Acts bound variant measurements.
/// @note The type IndexSourceLink is used for the measurements' source links.
template <typename detector_t, typename allocator_t>
inline auto createActsMeasurements(const detector_t& detector, const std::vector<traccc::measurement, allocator_t>& measurements){
    std::vector<Acts::BoundVariantMeasurement> measurementContainer;
    for (const traccc::measurement& m : measurements)
    {
        Acts::GeometryIdentifier moduleGeoId(detector.surface(m.surface_link).source);
        Index measurementIdx = measurementContainer.size();
        IndexSourceLink idxSourceLink{moduleGeoId, measurementIdx};
        measurementContainer.push_back(Acts::TracccPlugin::boundVariantMeasurement(m, Acts::SourceLink{idxSourceLink}));
    }
    return measurementContainer;
}

}
