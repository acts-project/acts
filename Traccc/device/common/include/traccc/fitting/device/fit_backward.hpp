/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/fitting/device/fit.hpp"
#include "traccc/fitting/status_codes.hpp"

namespace traccc::device {

template <typename fitter_t>
TRACCC_HOST_DEVICE inline void fit_backward(
    const global_index_t globalIndex, const typename fitter_t::config_type cfg,
    const fit_payload<fitter_t>& payload) {

    typename fitter_t::detector_type det(payload.det_data);

    vecmem::device_vector<const unsigned int> param_ids(payload.param_ids_view);
    vecmem::device_vector<unsigned int> param_liveness(
        payload.param_liveness_view);
    typename edm::track_container<
        typename fitter_t::detector_type::algebra_type>::device
        tracks(payload.tracks_view);

    if (globalIndex >= tracks.tracks.size()) {
        return;
    }

    const unsigned int param_id = param_ids.at(globalIndex);
    edm::track track = tracks.tracks.at(param_id);

    // Run fitting
    fitter_t fitter(det, payload.field_data, cfg);

    if (param_liveness.at(param_id) > 0u) {
        typename fitter_t::state fitter_state(
            track, tracks.states, tracks.measurements,
            *(payload.surfaces_view.ptr() + param_id),
            fitter.config().propagation, fitter.config().meas_calibration);

        kalman_fitter_status fit_status = fitter.smooth(fitter_state);

        fitter.update_statistics(fitter_state);

        // Assume that this branch is only called if the forward fit was
        // successfull (track param are alive)
        fitter.check_fitting_result(fitter_state, kalman_fitter_status::SUCCESS,
                                    fit_status);

        if (fit_status == kalman_fitter_status::SUCCESS) {
            track = fitter_state.m_fit_res;
        } else {
            param_liveness.at(param_id) = 0u;
        }

        // TODO: Grab the smoothed state for next it
    }
}

}  // namespace traccc::device
