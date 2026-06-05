/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/edm/device/sort_key.hpp"
#include "traccc/edm/track_parameters.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// (Event Data) Payload for the @c
/// traccc::device::fill_finding_propagation_sort_keys function
struct fill_finding_propagation_sort_keys_payload {
    /**
     * @brief View object to the vector of bound track parameters
     */
    bound_track_parameters_collection_types::const_view params_view;

    /**
     * @brief View object to the vector of boolean-like integers describing the
     * liveness of each parameter
     */
    vecmem::data::vector_view<const unsigned int> param_liveness_view;

    /**
     * @brief View object to the vector of sort keys
     */
    vecmem::data::vector_view<device::sort_key> keys_view;

    /**
     * @brief View object to the vector of parameter indices, which is the
     * output to the algorithm
     */
    vecmem::data::vector_view<unsigned int> ids_view;
};

/// Function used for fill key container
///
/// @param[in] globalIndex   The index of the current thread
/// @param[inout] payload      The function call payload
///
TRACCC_HOST_DEVICE inline void fill_finding_propagation_sort_keys(
    global_index_t globalIndex,
    const fill_finding_propagation_sort_keys_payload& payload);

}  // namespace traccc::device

// Include the implementation.
#include "./impl/fill_finding_propagation_sort_keys.ipp"
