/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/alpaka/clusterization/measurement_sorting_algorithm.hpp"

#include "../utils/get_queue.hpp"
#include "../utils/parallel_algorithms.hpp"
#include "../utils/thread_id.hpp"
#include "../utils/utils.hpp"

// Project include(s).
#include "traccc/clusterization/device/measurement_sorting.hpp"

// System include(s).
#include <memory_resource>

namespace traccc::alpaka {
namespace kernels {

/// Kernel wrapping @c traccc::device::fill_measurement_sort_keys
struct fill_measurement_sort_keys {
  /// @param[in]  acc               Alpaka accelerator object
  /// @param[in]  measurements_view View of the unsorted measurements
  /// @param[out] keys_view         View of the sorting keys
  /// @param[out] indices_view      View of the measurement indices
  ///
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc,
      const edm::measurement_collection::const_view measurements_view,
      vecmem::data::vector_view<device::measurement_sort_key_t> keys_view,
      vecmem::data::vector_view<unsigned int> indices_view) const {
    device::fill_measurement_sort_keys(
        details::thread_id1{acc}.getGlobalThreadId(), measurements_view,
        keys_view, indices_view);
  }
};  // struct fill_measurement_sort_keys

/// Kernel filling the output buffer with sorted measurements.
struct fill_sorted_measurements {
  /// @param[in] acc Alpaka accelerator object
  /// @param[in] input_view View of the input measurements
  /// @param[out] output_view View of the output measurements
  /// @param[in] sorted_indices_view View of the sorted measurement indices
  ///
  template <typename TAcc>
  ALPAKA_FN_ACC void operator()(
      TAcc const& acc, const edm::measurement_collection::const_view input_view,
      edm::measurement_collection::view output_view,
      const vecmem::data::vector_view<const unsigned int> sorted_indices_view)
      const {
    device::fill_sorted_measurements(
        details::thread_id1{acc}.getGlobalThreadId(), input_view, output_view,
        sorted_indices_view);
  }
};  // struct fill_sorted_measurements

}  // namespace kernels

measurement_sorting_algorithm::measurement_sorting_algorithm(
    const traccc::memory_resource& mr, const vecmem::copy& copy, queue& q,
    std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), m_mr{mr}, m_copy{copy}, m_queue{q} {}

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

  auto queue = details::get_queue(m_queue);

  // Sorting keys and index sequence.
  vecmem::data::vector_buffer<device::measurement_sort_key_t> keys(
      n_measurements, m_mr.main);
  vecmem::data::vector_buffer<unsigned int> indices(n_measurements, m_mr.main);
  m_copy.get().setup(keys)->wait();
  m_copy.get().setup(indices)->wait();

  static constexpr unsigned int BLOCK_SIZE = 256;
  const unsigned int n_blocks = (n_measurements + BLOCK_SIZE - 1) / BLOCK_SIZE;
  auto workDiv = makeWorkDiv<Acc>(n_blocks, BLOCK_SIZE);

  // Sort the indices by the sorting keys, with a radix sort.
  ::alpaka::exec<Acc>(queue, workDiv, kernels::fill_measurement_sort_keys{},
                      measurements_view, vecmem::get_data(keys),
                      vecmem::get_data(indices));
  details::sort_by_key(queue, m_mr, keys.ptr(), keys.ptr() + n_measurements,
                       indices.ptr());

  // Fill the output with the sorted measurements.
  ::alpaka::exec<Acc>(queue, workDiv, kernels::fill_sorted_measurements{},
                      measurements_view, vecmem::get_data(result),
                      vecmem::get_data(indices));

  // Return the sorted buffer.
  return result;
}

}  // namespace traccc::alpaka
