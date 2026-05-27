/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/io/csv/dfe.hpp"

// System include(s).
#include <cstdint>

namespace traccc::io::csv {

/// Type used in reading CSV data into memory
struct cell {

    uint64_t geometry_id = 0;
    uint64_t measurement_id = 0;
    uint32_t channel0 = 0;
    uint32_t channel1 = 0;
    float timestamp = 0.f;
    float value = 0.f;

    auto operator<=>(const cell& other) const = default;

    // geometry_id,measurement_id,channel0,channel1,timestamp,value
    DFE_NAMEDTUPLE(cell, geometry_id, measurement_id, channel0, channel1,
                   timestamp, value);
};

}  // namespace traccc::io::csv
