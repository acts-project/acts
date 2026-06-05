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

/// (Event Data) Payload for the @c
/// traccc::device::sort_updated_tracks function
struct sort_updated_tracks_payload {

    /**
     * @brief View object to the vector of relative number of shared
     * measurements
     */
    vecmem::data::vector_view<const traccc::scalar> rel_shared_view;

    /**
     * @brief View object to the vector of relative number of pvalues
     */
    vecmem::data::vector_view<const traccc::scalar> pvals_view;

    /**
     * @brief Whether to terminate the calculation
     */
    int* terminate;

    /**
     * @brief The number of updated tracks
     */
    unsigned int* n_updated_tracks;

    /**
     * @brief View object to the updated track
     */
    vecmem::data::vector_view<unsigned int> updated_tracks_view;
};

}  // namespace traccc::device
