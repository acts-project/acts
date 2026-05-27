/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/thread_id.hpp"
#include "reify_cluster_data.cuh"

// Project include(s).
#include "traccc/clusterization/device/reify_cluster_data.hpp"

namespace traccc::cuda::kernels {
__global__ void reify_cluster_data(
    vecmem::data::vector_view<const unsigned int> disjoint_set_view,
    traccc::edm::silicon_cluster_collection::view cluster_view) {

    device::reify_cluster_data(details::thread_id1{}.getGlobalThreadId(),
                               disjoint_set_view, cluster_view);
}
}  // namespace traccc::cuda::kernels
