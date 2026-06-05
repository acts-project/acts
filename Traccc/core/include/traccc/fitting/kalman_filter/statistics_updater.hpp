/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_state_collection.hpp"

namespace traccc {

/// Type unrolling functor to update the fitting quality
template <typename algebra_t>
struct statistics_updater {

    using scalar_type = detray::dscalar<algebra_t>;

    /// Update track fitting qualities (NDoF and Chi2)
    ///
    /// @param mask_group mask group that contains the mask of the current
    /// surface
    /// @param index mask index of the current surface
    /// @param fit_res fitting information such as NDoF or Chi2
    /// @param trk_state track state of the current surface
    TRACCC_HOST_DEVICE inline void operator()(
        typename edm::track_collection<algebra_t>::device::proxy_type& fit_res,
        const typename edm::track_state_collection<
            algebra_t>::const_device::const_proxy_type& trk_state,
        const edm::measurement_collection::const_device& measurements) {

        if (!trk_state.is_hole()) {

            // Measurement dimension
            const unsigned int D =
                measurements.at(trk_state.measurement_index()).dimensions();

            if (trk_state.is_smoothed()) {
                // NDoF = NDoF + number of coordinates per measurement
                fit_res.ndf() += static_cast<scalar_type>(D);
                fit_res.chi2() += trk_state.backward_chi2();
            }
        }
    }
};

}  // namespace traccc
