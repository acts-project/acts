/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "../sanity/contiguous_on.hpp"
#include "traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"

// Project include(s).
#include "traccc/edm/device/identity_projector.hpp"

// HIP include(s).
#include <hip/hip_runtime.h>

namespace traccc::hip {

bool combinatorial_kalman_filter_algorithm::input_is_valid(
    const edm::measurement_collection::const_view& measurements) const {
  static constexpr std::size_t GEOMID_INDEX = 6u;
  return is_contiguous_on<
      vecmem::device_vector<const detray::geometry::identifier>>(
      device::identity_projector{}, mr().main, copy(), stream(),
      measurements.template get<GEOMID_INDEX>());
}

}  // namespace traccc::hip
