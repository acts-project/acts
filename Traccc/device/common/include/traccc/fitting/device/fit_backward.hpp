/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"
#include "traccc/fitting/device/fit_payload.hpp"
#include "traccc/fitting/status_codes.hpp"

namespace traccc::device {

/// Function performing a backward fit iteration
template <typename fitter_t>
TRACCC_HOST_DEVICE inline void fit_backward(
    const global_index_t globalIndex, const typename fitter_t::config_type& cfg,
    const fit_payload& payload,
    const fit_tpayload<typename fitter_t::detector_type::const_view_type,
                       typename fitter_t::bfield_type,
                       typename fitter_t::surface_type>& tpayload) {
  typename fitter_t::detector_type det(tpayload.det);

  vecmem::device_vector<const unsigned int> param_ids(payload.track_indices);
  vecmem::device_vector<unsigned int> param_liveness(payload.track_liveness);
  typename edm::track_container<
      typename fitter_t::detector_type::algebra_type>::device
      tracks(payload.tracks);

  if (globalIndex >= tracks.tracks.size()) {
    return;
  }

  const unsigned int param_id =
      param_ids.size() == 0u ? globalIndex : param_ids.at(globalIndex);
  edm::track track = tracks.tracks.at(param_id);

  // Run fitting
  fitter_t fitter(det, tpayload.field, cfg);

  if ((param_liveness.empty() && track.constituent_links().size() > 0u) ||
      (param_id < param_liveness.size() && param_liveness.at(param_id) > 0u)) {
    typename fitter_t::state fitter_state(
        track, tracks.states, tracks.measurements,
        *(tpayload.surfaces.ptr() + param_id), fitter.config().propagation,
        fitter.config().meas_calibration);

    const kalman_fitter_status fit_status = fitter.smooth(fitter_state);

    fitter.update_statistics(fitter_state);

    // Assume that this branch is only called if the forward fit was
    // successfull (track param are alive)
    fitter.check_fitting_result(fitter_state, kalman_fitter_status::SUCCESS,
                                fit_status);

    if (!param_liveness.empty() &&
        fit_status != kalman_fitter_status::SUCCESS) {
      param_liveness.at(param_id) = 0u;
    }

    // TODO: Grab the smoothed state for next it
  }
}
}  // namespace traccc::device
