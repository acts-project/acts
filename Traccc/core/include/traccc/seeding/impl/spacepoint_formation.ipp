/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/measurement_helpers.hpp"

// Detray include(s).
#include <detray/geometry/tracking_surface.hpp>

namespace traccc::details {

template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline bool is_valid_measurement(
    const edm::measurement<measurement_backend_t>& meas) {
    // We use 2D (pixel) measurements only for spacepoint creation
    return (meas.dimensions() == 2u);
}

template <typename spacepoint_backend_t, typename detector_t,
          typename measurement_backend_t>
TRACCC_HOST_DEVICE inline void fill_pixel_spacepoint(
    edm::spacepoint<spacepoint_backend_t>& sp, const detector_t& det,
    const edm::measurement<measurement_backend_t>& meas,
    const typename detector_t::geometry_context gctx) {

    // Get the global position of this silicon pixel measurement.
    const detray::tracking_surface sf{det, meas.surface_link()};
    const detray::dpoint3D<typename detector_t::algebra_type> global =
        sf.local_to_global(
            gctx,
            edm::get_measurement_local<typename detector_t::algebra_type>(meas),
            {});

    // Fill the spacepoint with the global position and the measurement.
    sp.x() = static_cast<float>(getter::element(global, 0u));
    sp.y() = static_cast<float>(getter::element(global, 1u));
    sp.z() = static_cast<float>(getter::element(global, 2u));
    sp.radius_variance() = 0.f;
    sp.z_variance() = 0.f;
}

}  // namespace traccc::details
