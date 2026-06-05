/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../../utils/global_index.hpp"
#include "count_shared_measurements.cuh"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/jagged_device_vector.hpp>

// Thrust include(s).
#include <thrust/binary_search.h>
#include <thrust/execution_policy.h>

namespace traccc::cuda::kernels {

__global__ void count_shared_measurements(
    device::count_shared_measurements_payload payload) {

    vecmem::device_vector<const unsigned int> accepted_ids(
        payload.accepted_ids_view);

    auto globalIndex = threadIdx.x + blockIdx.x * blockDim.x;

    if (globalIndex >= accepted_ids.size()) {
        return;
    }

    vecmem::jagged_device_vector<const measurement_id_type> meas_ids(
        payload.meas_ids_view);
    vecmem::device_vector<const unsigned int> meas_id_to_unique_id(
        payload.meas_id_to_unique_id_view);
    vecmem::device_vector<const unsigned int> n_accepted_tracks_per_measurement(
        payload.n_accepted_tracks_per_measurement_view);
    vecmem::device_vector<unsigned int> n_shared(payload.n_shared_view);

    const unsigned int id = accepted_ids.at(globalIndex);

    for (const auto& meas_id : meas_ids[id]) {

        const auto unique_meas_idx = meas_id_to_unique_id.at(meas_id);

        if (n_accepted_tracks_per_measurement.at(unique_meas_idx) > 1) {
            vecmem::device_atomic_ref<unsigned int>(n_shared.at(id))
                .fetch_add(1u);
        }
    }
}
}  // namespace traccc::cuda::kernels
