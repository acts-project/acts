/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/clusterization/measurement_sorting_algorithm.hpp"

// Project include(s).
#include "traccc/clusterization/device/measurement_sorting.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/utils/copy.hpp>

// Thrust include(s).
#include <thrust/execution_policy.h>
#include <thrust/sort.h>

// System include(s).
#include <memory_resource>

namespace traccc::cuda {
namespace kernels {

/// Kernel wrapping @c traccc::device::fill_measurement_sort_keys
__global__ void fill_measurement_sort_keys(
    const edm::measurement_collection::const_view measurements_view,
    vecmem::data::vector_view<device::measurement_sort_key_t> keys_view,
    vecmem::data::vector_view<unsigned int> indices_view) {
  device::fill_measurement_sort_keys(
      details::global_index1(), measurements_view, keys_view, indices_view);
}

/// Kernel wrapping @c traccc::device::fill_sorted_measurements
__global__ void fill_sorted_measurements(
    const edm::measurement_collection::const_view input_view,
    edm::measurement_collection::view output_view,
    const vecmem::data::vector_view<const unsigned int> sorted_indices_view) {
  device::fill_sorted_measurements(details::global_index1(), input_view,
                                   output_view, sorted_indices_view);
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

  // Get the number of measurements.
  edm::measurement_collection::const_view::size_type n_measurements = 0u;
  if (m_mr.host) {
    const vecmem::async_size size =
        m_copy.get().get_size(measurements_view, *(m_mr.host));
    n_measurements = size.get();
  } else {
    n_measurements = m_copy.get().get_size(measurements_view);
  }

  // Create the output buffer.
  output_type result{measurements_view.capacity(), m_mr.main,
                     vecmem::data::buffer_type::resizable};
  m_copy.get().setup(result)->ignore();
  if (n_measurements == 0) {
    return result;
  }
  m_copy.get()(measurements_view.size(), result.size())->ignore();

  // Get a convenience variable for the stream that we'll be using.
  cudaStream_t stream = details::get_stream(m_stream);
  // Set up the Thrust execution policy.
  auto policy =
      thrust::cuda::par_nosync(std::pmr::polymorphic_allocator(&(m_mr.main)))
          .on(stream);

  // Sorting keys and index sequence.
  vecmem::data::vector_buffer<device::measurement_sort_key_t> keys(
      n_measurements, m_mr.main);
  vecmem::data::vector_buffer<unsigned int> indices(n_measurements, m_mr.main);
  m_copy.get().setup(keys)->ignore();
  m_copy.get().setup(indices)->ignore();

  static constexpr unsigned int BLOCK_SIZE = 256;
  const unsigned int n_blocks = (n_measurements + BLOCK_SIZE - 1) / BLOCK_SIZE;

  // Sort the indices by the sorting keys, with a radix sort.
  kernels::fill_measurement_sort_keys<<<n_blocks, BLOCK_SIZE, 0, stream>>>(
      measurements_view, keys, indices);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
  thrust::sort_by_key(policy, keys.ptr(), keys.ptr() + n_measurements,
                      indices.ptr());

  // Fill the output with the sorted measurements.
  kernels::fill_sorted_measurements<<<n_blocks, BLOCK_SIZE, 0, stream>>>(
      measurements_view, result, indices);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());

  // Return the sorted buffer.
  return result;
}

}  // namespace traccc::cuda
