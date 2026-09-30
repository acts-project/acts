/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/seeding/device/silicon_pixel_spacepoint_formation_algorithm.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>

// System include(s).
#include <vector>

namespace traccc::device {

silicon_pixel_spacepoint_formation_algorithm::
    silicon_pixel_spacepoint_formation_algorithm(
        const traccc::memory_resource& mr, const vecmem::copy& copy,
        std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), algorithm_base(mr, copy) {}

auto silicon_pixel_spacepoint_formation_algorithm::operator()(
    const detector_buffer& det,
    const edm::measurement_collection::const_view& measurements) const
    -> output_type {
  // Get the number of measurements. In an asynchronous way if possible.
  edm::measurement_collection::const_view::size_type n_measurements = 0u;
  if (mr().host) {
    vecmem::async_size size = copy().get_size(measurements, *(mr().host));
    // Block or suspend execution until the size is available.
    await(size);
    n_measurements = size.unsafe_get();
  } else {
    n_measurements = copy().get_size(measurements);
  }

  // If there are no measurements, return right away.
  if (n_measurements == 0) {
    return {};
  }

  // Flag every measurement that a spacepoint is made out of, and turn those
  // flags into a prefix sum. The prefix sum gives every measurement the index
  // of its spacepoint, which makes the output independent of the thread order.
  vecmem::data::vector_buffer<unsigned int> spacepoint_flags(n_measurements,
                                                             mr().main);
  copy().setup(spacepoint_flags)->ignore();
  vecmem::data::vector_view<unsigned int> spacepoint_flags_view(
      spacepoint_flags);
  count_spacepoints_kernel(
      {n_measurements, measurements, spacepoint_flags_view});
  scan_spacepoint_flags(spacepoint_flags_view);

  // The scan is inclusive, so its last element is the number of spacepoints.
  // We need that number on the host to size the output buffer, so we have to
  // wait for the scan here.
  assert(mr().host != nullptr);
  vecmem::vector<unsigned int> n_spacepoints_host(mr().host);
  // TODO: Yield during this asynchronous copy.
  copy()(vecmem::data::vector_view<const unsigned int>(
             1u, spacepoint_flags_view.ptr() + n_measurements - 1u),
         n_spacepoints_host)
      ->wait();
  const unsigned int n_spacepoints = n_spacepoints_host.at(0);

  // If no measurement produces a spacepoint, return right away.
  if (n_spacepoints == 0) {
    return {};
  }

  // Create the result buffer.
  edm::spacepoint_collection::buffer spacepoints(n_spacepoints, mr().main);
  copy().setup(spacepoints)->ignore();

  // Launch the spacepoint formation kernel.
  const vecmem::data::vector_view<const unsigned int> spacepoint_index_view(
      spacepoint_flags_view);
  form_spacepoints_kernel(
      {n_measurements, det, measurements, spacepoint_index_view, spacepoints});

  // Return the reconstructed spacepoints.
  return spacepoints;
}

}  // namespace traccc::device
