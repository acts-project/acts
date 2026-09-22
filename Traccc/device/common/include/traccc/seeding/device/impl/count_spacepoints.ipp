/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/seeding/detail/spacepoint_formation.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device {

TRACCC_HOST_DEVICE inline void count_spacepoints(
    const global_index_t globalIndex,
    const edm::measurement_collection::const_view& measurements_view,
    vecmem::data::vector_view<unsigned int> spacepoint_flags_view) {
  // Set up the input container(s).
  const edm::measurement_collection::const_device measurements(
      measurements_view);

  // Check if anything needs to be done
  if (globalIndex >= measurements.size()) {
    return;
  }

  // Set up the output container(s).
  vecmem::device_vector<unsigned int> spacepoint_flags(spacepoint_flags_view);

  // Flag the measurement if a spacepoint is created out of it.
  spacepoint_flags.at(globalIndex) =
      (traccc::details::is_valid_measurement(measurements.at(globalIndex))
           ? 1u
           : 0u);
}

}  // namespace traccc::device
