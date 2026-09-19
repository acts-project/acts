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
#include "traccc/geometry/detector.hpp"
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
  device::estimate_track_params(details::global_index1(), config, measurements,
                                spacepoints, seeds, bfield, params_view);
}

/// CUDA kernel for running @c traccc::device::estimate_track_params with a
/// detector
template <typename detector_t, typename bfield_t>
__global__ void estimate_track_params_with_detector(
    const track_params_estimation_config config,
    typename detector_t::view det_view,
    edm::measurement_collection::const_view measurements,
    edm::spacepoint_collection::const_view spacepoints,
    edm::seed_collection::const_view seeds, const bfield_t bfield,
    bound_track_parameters_collection_types::view params_view)
  requires(traccc::is_detector_traits<detector_t>)
{
  const typename detector_t::device det{det_view};
  device::estimate_track_params(details::global_index1(), config, det,
                                measurements, spacepoints, seeds, bfield,
                                params_view);
}

}  // namespace kernels

seed_parameter_estimation_algorithm::seed_parameter_estimation_algorithm(
    const track_params_estimation_config& config,
    const traccc::memory_resource& mr, const vecmem::copy& copy,
    const stream_wrapper& str, std::unique_ptr<const Logger> logger,
    await_function_type await_func)
    : device::seed_parameter_estimation_algorithm(config, mr, copy,
                                                  std::move(logger)),
      cuda::algorithm_base(str, std::move(await_func)) {}

void seed_parameter_estimation_algorithm::estimate_seed_params_kernel(
    const struct estimate_seed_params_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 4;
  const unsigned int n_blocks = (payload.n_seeds + n_threads - 1) / n_threads;
  magnetic_field_visitor<bfield_type_list<scalar>>(
      payload.bfield, [&]<typename bfield_view_t>(const bfield_view_t& bfield) {
        if (payload.detector == nullptr) {
          kernels::estimate_track_params<<<n_blocks, n_threads, 0,
                                           details::get_stream(stream())>>>(
              payload.config, payload.measurements, payload.spacepoints,
              payload.seeds, bfield, payload.params);
        } else {
          detector_buffer_visitor<detector_type_list>(
              *(payload.detector),
              [&]<typename detector_traits_t>(
                  const typename detector_traits_t::view& det) {
                kernels::estimate_track_params_with_detector<detector_traits_t>
                    <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(
                        payload.config, det, payload.measurements,
                        payload.spacepoints, payload.seeds, bfield,
                        payload.params);
              });
        }
      });
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

}  // namespace traccc::cuda
