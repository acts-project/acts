/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"
#include "traccc/fitting/status_codes.hpp"

namespace traccc::host {

/// Run the Kalman smoother over the collected surfaces
template <typename fitter_t>
TRACCC_HOST inline void kalman_smoother(
    const typename fitter_t::config_type cfg,
    const typename fitter_t::detector_type& det,
    const typename fitter_t::bfield_type& field,
    typename edm::track_container<
        typename fitter_t::detector_type::algebra_type>::host& track_container,
    vecmem::data::jagged_vector_view<
        typename fitter_t::detector_type::surface_type>
        surfaces_view) {

    using algebra_t = typename fitter_t::detector_type::algebra_type;

    // Run fitting
    fitter_t fitter(det, field, cfg);

    // Create the input container(s).
    typename edm::track_container<algebra_t>::data tracks_data{track_container};
    typename edm::track_container<algebra_t>::view tracks_view{tracks_data};
    typename edm::track_container<algebra_t>::device device_track_container{
        tracks_view};
    edm::measurement_collection::const_device measurements{
        track_container.measurements};

    for (unsigned int trk_idx = 0u; trk_idx < track_container.tracks.size();
         trk_idx++) {

        // Device(!) container proxy type
        auto track = device_track_container.tracks.at(trk_idx);
        if (track.constituent_links().size() == 0u) {
            continue;
        }
        TRACCC_INFO_HOST("Fitting track " << trk_idx << " ("
                                          << track.constituent_links().size()
                                          << " track states)");
        TRACCC_VERBOSE_HOST("" << track_container.tracks.at(trk_idx));

        vecmem::data::vector_view<
            typename fitter_t::detector_type::surface_type>
            sf_view{*(surfaces_view.ptr() + trk_idx)};

        typename fitter_t::state fitter_state(
            track, device_track_container.states, measurements, sf_view,
            fitter.config().propagation, fitter.config().meas_calibration);

        const kalman_fitter_status fit_status = fitter.smooth(fitter_state);

        fitter.update_statistics(fitter_state);

        // Assume that this branch is only called if the forward fit was
        // successfull (track param are alive)
        fitter.check_fitting_result(fitter_state, kalman_fitter_status::SUCCESS,
                                    fit_status);
    }
}

}  // namespace traccc::host
