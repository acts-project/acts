/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once
// Project includes
#include "traccc/definitions/primitives.hpp"

// Acts includes
#include <Acts/Geometry/GeometryHierarchyMap.hpp>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace traccc {

/// Type describing the conditions configuration of a detector module
struct conditions_data_config {
    vector2 shift{0.f, 0.f};
};

using conditions_config = Acts::GeometryHierarchyMap<conditions_data_config>;

}  // namespace traccc
