/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::edm {

template <typename algebra_t>
TRACCC_HOST_DEVICE
    typename track_state_collection<algebra_t>::device::object_type
    make_track_state(const measurement_collection::const_device& measurements,
                     unsigned int mindex) {

    // Create the result object.
    typename track_state_collection<algebra_t>::device::object_type state;

    // Set it not to be a hole by default, with the appropriate (measurement)
    // index.
    state.set_hole(false);
    state.measurement_index() = mindex;

    // Set the correct surface link for the track parameters.
    state.filtered_params().set_surface_link(
        measurements.at(mindex).surface_link());
    state.smoothed_params().set_surface_link(
        measurements.at(mindex).surface_link());

    // Return the initialized state.
    return state;
}

}  // namespace traccc::edm
