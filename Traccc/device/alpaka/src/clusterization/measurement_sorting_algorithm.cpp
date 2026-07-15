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

// Project include(s).
#include "traccc/clusterization/device/geo_id_based_sorter.hpp"
#include "traccc/clusterization/device/sorting_index_filler.hpp"

// System include(s).
#include <memory_resource>

namespace traccc::alpaka {
namespace kernels {

/// Kernel filling the output buffer with sorted measurements.
struct fill_sorted_measurements {
    /// @param[in] acc Alpaka accelerator object
    /// @param[in] input_view View of the input measurements
    /// @param[out] output_view View of the output measurements
    /// @param[in] sorted_indices_view View of the sorted measurement indices
    ///
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        const edm::measurement_collection::const_view input_view,
        edm::measurement_collection::view output_view,
        const vecmem::data::vector_view<const unsigned int> sorted_indices_view)
        const {

        // Create the device objects.
        const edm::measurement_collection::const_device input{input_view};
        edm::measurement_collection::device output{output_view};
        const vecmem::device_vector<const unsigned int> sorted_indices{
            sorted_indices_view};

        // Stop early if we can.
        const unsigned int index = details::thread_id1{acc}.getGlobalThreadId();
        if (index >= input.size()) {
            return;
        }

        // Copy one measurement into the correct position.
        output.at(index) = input.at(sorted_indices.at(index));
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

    // Get a convenience variable for the queue that we'll be using.
    auto queue = details::get_queue(m_queue);

    // Create a device container on top of the view.
    const edm::measurement_collection::const_device measurements{
        measurements_view};

    // Create a vector of measurement indices, which would be sorted.
    vecmem::data::vector_buffer<unsigned int> indices(
        measurements_view.capacity(), m_mr.main);
    m_copy.get().setup(indices)->wait();
    details::for_each(queue, m_mr, indices.ptr(),
                      indices.ptr() + indices.capacity(),
                      device::sorting_index_filler{indices});

    // Sort the indices according to the surface identifiers of the
    // measurements.
    details::sort(queue, m_mr, indices.ptr(),
                  indices.ptr() + indices.capacity(),
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
    auto workDiv = makeWorkDiv<Acc>(n_blocks, BLOCK_SIZE);
    ::alpaka::exec<Acc>(queue, workDiv, kernels::fill_sorted_measurements{},
                        measurements_view, vecmem::get_data(result),
                        vecmem::get_data(indices));

    // Return the sorted buffer.
    return result;
}

}  // namespace traccc::alpaka
