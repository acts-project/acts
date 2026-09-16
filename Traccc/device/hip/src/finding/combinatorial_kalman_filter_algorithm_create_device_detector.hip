/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../utils/hip_error_handling.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

// System include(s).
#include <new>
#include <type_traits>
#include <utility>

namespace traccc::hip {
namespace kernels {

/// HIP kernel constructing the device detector in place, in global memory
///
/// The detector is built once rather than in every thread, so that the
/// propagation kernels only carry a pointer to it.
///
template <typename detector_t>
__global__ void create_device_detector(typename detector_t::const_view_type in,
                                       detector_t* out) {
  unsigned int thread_id = blockIdx.x * blockDim.x + threadIdx.x;

  if (thread_id == 0) {
    new (out) detector_t(in);
  }
}

}  // namespace kernels

move_only_any combinatorial_kalman_filter_algorithm::create_device_detector(
    const detector_buffer& det) const {
  return detector_buffer_visitor<detector_type_list>(
      det, [&]<typename detector_traits_t>(
               const typename detector_traits_t::view& det_view) {
        using device_detector_t = typename detector_traits_t::device;

        static_assert(
            std::is_trivially_destructible_v<device_detector_t>,
            "the detector is placement-constructed and never destroyed, so it "
            "must not need a destructor");

        vecmem::data::vector_buffer<device_detector_t> device_detector_buffer(
            1, mr().main);
        hipLaunchKernelGGL(kernels::create_device_detector<device_detector_t>,
                           1, 1, 0, details::get_stream(stream()), det_view,
                           device_detector_buffer.ptr());
        TRACCC_HIP_ERROR_CHECK(hipGetLastError());
        return move_only_any{std::move(device_detector_buffer)};
      });
}

}  // namespace traccc::hip
