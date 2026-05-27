/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::device {

TRACCC_HOST_DEVICE inline void reify_cluster_data(
    global_index_t thread_id,
    vecmem::data::vector_view<const unsigned int> disjoint_set_view,
    traccc::edm::silicon_cluster_collection::view cluster_view) {

    // Create the device objects.
    const vecmem::device_vector<const unsigned int> disjoint_set(
        disjoint_set_view);
    traccc::edm::silicon_cluster_collection::device clusters(cluster_view);

    // Fill the output container.
    if (thread_id < disjoint_set.size()) {
        clusters.cell_indices()
            .at(disjoint_set.at(thread_id))
            .push_back(thread_id);
    }
}

}  // namespace traccc::device
