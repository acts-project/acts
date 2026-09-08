/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "../utils/hip_error_handling.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"
#include "traccc/hip/clusterization/measurement_sorting_algorithm.hpp"

// Project include(s).
#include "traccc/clusterization/device/geo_id_based_sorter.hpp"
#include "traccc/clusterization/device/sorting_index_filler.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>

// (Roc)Thrust include(s).
#include <thrust/for_each.h>
#include <thrust/sort.h>

// HIP include(s).
#include <hip/hip_runtime.h>

// System include(s).
#include <memory_resource>

namespace traccc::hip {
namespace kernels {

/// Kernel filling the output buffer with sorted measurements.
///
/// @param[in] input_view View of the input measurements
/// @param[out] output_view View of the output measurements
/// @param[in] sorted_indices_view View of the sorted measurement indices
///
__global__ void fill_sorted_measurements(
    const edm::measurement_collection::const_view input_view,
    edm::measurement_collection::view output_view,
    const vecmem::data::vector_view<const unsigned int> sorted_indices_view) {
  // Create the device objects.
  const edm::measurement_collection::const_device input{input_view};
  edm::measurement_collection::device output{output_view};
  const vecmem::device_vector<const unsigned int> sorted_indices{
      sorted_indices_view};

  // Stop early if we can.
  const unsigned int index = details::thread_id1{}.getGlobalThreadId();
  if (index >= input.size()) {
    return;
  }

  // Copy one measurement into the correct position.
  output.at(index) = input.at(sorted_indices.at(index));
}

}  // namespace kernels

measurement_sorting_algorithm::measurement_sorting_algorithm(
    const traccc::memory_resource& mr, const vecmem::copy& copy,
    const stream_wrapper& str, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), m_mr{mr}, m_copy{copy}, m_stream{str} {}

measurement_sorting_algorithm::output_type
measurement_sorting_algorithm::operator()(
    const edm::measurement_collection::const_view& measurements_view) const {
  // Exit early if there are no measurements.
  if (measurements_view.capacity() == 0) {
    return {};
  }

  // Get a convenience variable for the stream that we'll be using.
  hipStream_t stream = details::get_stream(m_stream);
  // Set up the Thrust execution policy. Unfortunately this has to be done in
  // a backend specify way. :-(
#if defined(__HIP_PLATFORM_AMD__)
  auto policy =
      thrust::hip::par_nosync(std::pmr::polymorphic_allocator(&(m_mr.main)))
          .on(stream);
#elif defined(__HIP_PLATFORM_NVIDIA__)
  auto policy =
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(m_mr.main)))
          .on(stream);
#else
#error "Unrecognised HIP platform."
#endif

  // Create a device container on top of the view.
  const edm::measurement_collection::const_device measurements{
      measurements_view};

  // Create a vector of measurement indices, which would be sorted.
  vecmem::data::vector_buffer<unsigned int> indices(
      measurements_view.capacity(), m_mr.main);
  m_copy.get().setup(indices)->ignore();
  thrust::for_each(policy, indices.ptr(), indices.ptr() + indices.capacity(),
                   device::sorting_index_filler{indices});

  // Sort the indices according to the surface identifiers of the
  // measurements.
  thrust::sort(policy, indices.ptr(), indices.ptr() + indices.capacity(),
               device::geo_id_based_sorter{measurements.surface_link()});

  // Create the output buffer.
  output_type result{measurements_view.capacity(), m_mr.main,
                     vecmem::data::buffer_type::resizable};
  m_copy.get().setup(result)->ignore();
  m_copy.get()(measurements_view.size(), result.size())->ignore();

  // Fill it with the sorted measurements.
  static constexpr unsigned int BLOCK_SIZE = 256;
  const unsigned int n_blocks =
      (measurements_view.capacity() + BLOCK_SIZE - 1) / BLOCK_SIZE;
  hipLaunchKernelGGL(kernels::fill_sorted_measurements, n_blocks, BLOCK_SIZE, 0,
                     stream, measurements_view, result, indices);
  TRACCC_HIP_ERROR_CHECK(hipGetLastError());

  // Return the sorted buffer.
  return result;
}

}  // namespace traccc::hip
