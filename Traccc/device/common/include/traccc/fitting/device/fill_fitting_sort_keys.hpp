/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"
#include "traccc/edm/device/sort_key.hpp"

// Project include(s).
#include "traccc/edm/track_collection.hpp"

namespace traccc::device {

/// Function used to fill key container
///
/// @param[in] globalIndex   The index of the current thread
/// @param[in] track_candidates_view The input track candidates
/// @param[out] keys_view    The key values
/// @param[out] ids_view     The param ids
///
TRACCC_HOST_DEVICE inline void fill_fitting_sort_keys(
    global_index_t globalIndex,
    const edm::track_collection<default_algebra>::const_view&
        track_candidates_view,
    vecmem::data::vector_view<device::sort_key> keys_view,
    vecmem::data::vector_view<unsigned int> ids_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/fitting/device/impl/fill_fitting_sort_keys.ipp"
