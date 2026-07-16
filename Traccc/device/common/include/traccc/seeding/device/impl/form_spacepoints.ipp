/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/seeding/detail/spacepoint_formation.hpp"

// System include(s).
#include <cassert>

namespace traccc::device {

template <typename detector_t>
TRACCC_HOST_DEVICE inline void form_spacepoints(
    const global_index_t globalIndex, typename detector_t::view det_view,
    const edm::measurement_collection::const_view& measurements_view,
    edm::spacepoint_collection::view spacepoints_view) {

    // Set up the input container(s).
    const edm::measurement_collection::const_device measurements(
        measurements_view);

    // Check if anything needs to be done
    if (globalIndex >= measurements.size()) {
        return;
    }

    // Create the tracking geometry
    typename detector_t::device det(det_view);

    // Set up the output container(s).
    edm::spacepoint_collection::device spacepoints(spacepoints_view);

    const edm::measurement meas = measurements.at(globalIndex);

    // Fill the spacepoint using the common function.
    if (details::is_valid_measurement(meas)) {
        const edm::spacepoint_collection::device::size_type i =
            spacepoints.push_back_default();
        edm::spacepoint sp = spacepoints.at(i);
        traccc::details::fill_pixel_spacepoint(sp, det, meas);
        sp.measurement_index_1() = globalIndex;
        sp.measurement_index_2() =
            edm::spacepoint_collection::device::INVALID_MEASUREMENT_INDEX;
    }
}

}  // namespace traccc::device
