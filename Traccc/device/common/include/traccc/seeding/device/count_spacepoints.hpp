/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Function flagging the measurements that a spacepoint is made out of
///
/// The flags are the input of the prefix sum that gives every measurement its
/// spacepoint index. @c traccc::device::form_spacepoints uses that index, so
/// that the order of the spacepoints does not depend on the thread order.
///
/// @param[in] globalIndex           The index for the current thread
/// @param[in] measurements_view     Collection of measurements
/// @param[out] spacepoint_flags_view One flag for every measurement
///
TRACCC_HOST_DEVICE inline void count_spacepoints(
    global_index_t globalIndex,
    const edm::measurement_collection::const_view& measurements_view,
    vecmem::data::vector_view<unsigned int> spacepoint_flags_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/count_spacepoints.ipp"
