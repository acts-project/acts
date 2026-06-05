/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"

// VecMem include(s).
#include <vecmem/containers/data/jagged_vector_view.hpp>
#include <vecmem/containers/data/vector_view.hpp>

// System include(s).
#include <cstddef>

namespace traccc::device {

/// (Event Data) Payload for the @c traccc::device::sort_tracks_per_measurement
/// function
struct sort_tracks_per_measurement_payload {

    /**
     * @brief View object to the tracks per measurement
     */
    vecmem::data::jagged_vector_view<unsigned int> tracks_per_measurement_view;
};

}  // namespace traccc::device
