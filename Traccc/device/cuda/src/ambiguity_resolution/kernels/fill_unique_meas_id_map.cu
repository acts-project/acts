/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "fill_unique_meas_id_map.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

namespace traccc::cuda::kernels {

__global__ void fill_unique_meas_id_map(
    device::fill_unique_meas_id_map_payload payload) {

    vecmem::device_vector<const measurement_id_type> unique_meas(
        payload.unique_meas_view);

    const auto globalIndex = details::global_index1();
    if (globalIndex >= unique_meas.size()) {
        return;
    }

    vecmem::device_vector<unsigned int> meas_id_to_unique_id(
        payload.meas_id_to_unique_id_view);

    auto meas_id = unique_meas.at(globalIndex);
    meas_id_to_unique_id.at(meas_id) = globalIndex;
}

}  // namespace traccc::cuda::kernels
