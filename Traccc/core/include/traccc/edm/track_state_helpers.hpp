/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_state_collection.hpp"

namespace traccc::edm {

/// Create a track state with default values.
///
/// @param measurements The collection of measurements to use for initialization
/// @param mindex       The index of the measurement to associate with the state
///
/// @return A track state object initialized with default values
///
template <typename algebra_t>
TRACCC_HOST_DEVICE
    typename track_state_collection<algebra_t>::device::object_type
    make_track_state(const measurement_collection::const_device& measurements,
                     unsigned int mindex);

}  // namespace traccc::edm

// Include the implementation.
#include "traccc/edm/impl/track_state_helpers.ipp"
