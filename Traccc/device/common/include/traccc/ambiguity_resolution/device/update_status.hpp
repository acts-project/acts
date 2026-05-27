/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Event Data) Payload for the @c
/// traccc::device::update_status function
struct update_status_payload {

    /**
     * @brief Whether to terminate the calculation
     */
    int* terminate;

    /**
     * @brief The number of accepted tracks
     */
    unsigned int* n_accepted;

    /**
     * @brief The number of updated tracks
     */
    unsigned int* n_updated_tracks;

    /**
     * @brief View object to the inverted ids
     */
    vecmem::data::vector_view<const unsigned int> temp_sorted_ids_view;

    /**
     * @brief View object to the sorted track
     */
    vecmem::data::vector_view<unsigned int> sorted_ids_view;

    /**
     * @brief View object to the updated track
     */
    vecmem::data::vector_view<unsigned int> updated_tracks_view;

    /**
     * @brief View object to the whether track id is updated
     */
    vecmem::data::vector_view<int> is_updated_view;

    /**
     * @brief View object to the vector of number of shared measurements
     */
    vecmem::data::vector_view<const unsigned int> n_shared_view;

    /**
     * @brief The number of max shared
     */
    unsigned int* max_shared;
};

}  // namespace traccc::device
