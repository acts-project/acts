/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/magnetic_field_types.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/seeding/seed_parameter_estimation_algorithm.hpp"

// Project include(s).
#include "traccc/seeding/device/estimate_track_params.hpp"

namespace traccc::cuda {
namespace kernels {

/// CUDA kernel for running @c traccc::device::estimate_track_params
template <typename bfield_t>
__global__ void estimate_track_params(
    const track_params_estimation_config config,
    edm::measurement_collection::const_view measurements,
    edm::spacepoint_collection::const_view spacepoints,
    edm::seed_collection::const_view seeds, const bfield_t bfield,
    bound_track_parameters_collection_types::view params_view) {

    device::estimate_track_params(details::global_index1(), config,
                                  measurements, spacepoints, seeds, bfield,
                                  params_view);
}

}  // namespace kernels

seed_parameter_estimation_algorithm::seed_parameter_estimation_algorithm(
    const track_params_estimation_config& config,
    const traccc::memory_resource& mr, vecmem::copy& copy, cuda::stream& str,
    std::unique_ptr<const Logger> logger)
    : device::seed_parameter_estimation_algorithm(config, mr, copy,
                                                  std::move(logger)),
      cuda::algorithm_base(str) {}

void seed_parameter_estimation_algorithm::estimate_seed_params_kernel(
    const struct estimate_seed_params_kernel_payload& payload) const {

    const unsigned int n_threads = warp_size() * 4;
    const unsigned int n_blocks = (payload.n_seeds + n_threads - 1) / n_threads;
    magnetic_field_visitor<bfield_type_list<scalar>>(
        payload.bfield,
        [&]<typename bfield_view_t>(const bfield_view_t& bfield) {
            kernels::estimate_track_params<<<n_blocks, n_threads, 0,
                                             details::get_stream(stream())>>>(
                payload.config, payload.measurements, payload.spacepoints,
                payload.seeds, bfield, payload.params);
        });
    TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

}  // namespace traccc::cuda
