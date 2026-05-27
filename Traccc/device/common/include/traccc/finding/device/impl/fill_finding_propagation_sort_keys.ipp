/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::device {

TRACCC_HOST_DEVICE inline void fill_finding_propagation_sort_keys(
    const global_index_t globalIndex,
    const fill_finding_propagation_sort_keys_payload& payload) {

    const bound_track_parameters_collection_types::const_device params(
        payload.params_view);
    const vecmem::device_vector<const unsigned int> param_liveness(
        payload.param_liveness_view);

    // Keys
    vecmem::device_vector<device::sort_key> keys_device(payload.keys_view);

    // Param id
    vecmem::device_vector<unsigned int> ids_device(payload.ids_view);

    if (globalIndex >= keys_device.size()) {
        return;
    }

    /*
     * Adding a large constant factor to any dead tracks will ensure that they
     * all end up at the end of the array, and so they will produce minimal
     * thread divergence.
     */
    keys_device.at(globalIndex) =
        device::get_sort_key(params.at(globalIndex)) +
        (param_liveness.at(globalIndex) == 0u ? 10000.0f : 0);
    ids_device.at(globalIndex) = globalIndex;
}

}  // namespace traccc::device
