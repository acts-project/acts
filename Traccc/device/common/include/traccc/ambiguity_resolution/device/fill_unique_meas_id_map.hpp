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
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Event Data) Payload for the @c traccc::device::fill_unique_meas_id_map
/// function
struct fill_unique_meas_id_map_payload {

    /**
     * @brief View object to the unique measurement ids
     */
    vecmem::data::vector_view<const measurement_id_type> unique_meas_view;

    /**
     * @brief View object to the meas id to unique id map
     */
    vecmem::data::vector_view<unsigned int> meas_id_to_unique_id_view;
};

}  // namespace traccc::device
