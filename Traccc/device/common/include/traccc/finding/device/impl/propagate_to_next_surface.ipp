/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// @TODO: Remove this once alpaka compiles without warning.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic warning "-Wmaybe-uninitialized"
#endif

// Project include(s).
#include "traccc/finding/details/combinatorial_kalman_filter_types.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"

// Detray include(s).
#include <detray/utils/tuple_helpers.hpp>

namespace traccc::device {

template <typename propagator_t, typename bfield_t>
TRACCC_HOST_DEVICE inline void propagate_to_next_surface(
    const global_index_t globalIndex, const finding_config& cfg,
    const propagate_to_next_surface_payload<propagator_t, bfield_t>& payload) {

    using algebra_t = typename propagator_t::detector_type::algebra_type;
    using scalar_t = detray::dscalar<algebra_t>;

    if (globalIndex >= payload.n_in_params) {
        return;
    }

    // Theta id
    vecmem::device_vector<const unsigned int> param_ids(payload.param_ids_view);

    const unsigned int param_id = param_ids.at(globalIndex);

    // Links
    vecmem::device_vector<const candidate_link> links(payload.links_view);

    const unsigned int link_idx = payload.prev_links_idx + param_id;
    const candidate_link link = links.at(link_idx);
    assert(link.step == payload.step);
    const unsigned int n_cands = link.step + 1 - link.n_skipped;

    // Parameter liveness
    vecmem::device_vector<unsigned int> params_liveness(
        payload.params_liveness_view);

    // tips
    vecmem::device_vector<unsigned int> tips(payload.tips_view);
    vecmem::device_vector<unsigned int> tip_lengths(payload.tip_lengths_view);

    // Detector
    typename propagator_t::detector_type det(payload.det_data);

    // Parameters
    bound_track_parameters_collection_types::device params(payload.params_view);

    if (params_liveness.at(param_id) == 0u) {
        return;
    }

    // Input bound track parameter
    const bound_track_parameters<> in_par = params.at(param_id);

    // Create propagator
    detray::propagation::config prop_cfg{cfg.propagation};
    prop_cfg.navigation.estimate_scattering_noise = false;
    propagator_t propagator(prop_cfg);

    // Create propagator state
    typename propagator_t::state propagation(in_par, payload.field_data, det);
    propagation.set_particle(
        detail::correct_particle_hypothesis(cfg.ptc_hypothesis, in_par));
    propagation.stepping()
        .template set_constraint<detray::step::constraint::e_accuracy>(
            cfg.propagation.stepping.step_constraint);

    // Actor states
    // Pathlimit aborter
    typename detray::actor::pathlimit_aborter<scalar_t>::state aborter_state{};
    // Parameter updater
    typename detray::actor::parameter_updater_state<algebra_t> updater_state{
        prop_cfg, in_par};
    // CKF-interactor
    traccc::details::ckf_interactor_t::state interactor_state;
    // Momentum aborter
    typename detray::actor::momentum_aborter<scalar_t>::state
        momentum_aborter_state{};
    // CKF aborter
    typename ckf_aborter::state ckf_aborter_state{};

    /*
     * If we are running the MBF smoother, we need to accumulate the Jacobians
     * between the two sensitives multiplicatively. To this end, we ask the
     * parameter transporter to multiply the Jacobians into this matrix, which
     * is set to the multiplicative identity.
     */
    if (cfg.run_mbf_smoother) {
        assert(payload.tmp_jacobian_ptr != nullptr);

        payload.tmp_jacobian_ptr[param_id] =
            matrix::identity<bound_matrix<algebra_t>>();
        updater_state.set_full_jacobian(&payload.tmp_jacobian_ptr[param_id]);
    }

    // Notify the KF and material interaction only at the first propagation
    // initialization. For all subsequent initializations, they will have run
    // already on the previous step
    updater_state.notify_on_initial(link.step == 0);
    momentum_aborter_state.min_pT(static_cast<scalar_t>(cfg.min_pT));
    momentum_aborter_state.min_p(static_cast<scalar_t>(cfg.min_p));
    ckf_aborter_state.min_step_length = cfg.min_step_length_for_next_surface;
    ckf_aborter_state.max_count = cfg.max_step_counts_for_next_surface;

    // Propagate to the next surface
    propagator.propagate(
        propagation, detray::tie(aborter_state, updater_state, interactor_state,
                                 momentum_aborter_state, ckf_aborter_state));

    // If a surface found, add the parameter for the next step
    if (ckf_aborter_state.success) {
        assert(propagation.navigation().is_on_sensitive());
        assert(!updater_state.bound_params().is_invalid());

        params[param_id] = updater_state.bound_params();
        params_liveness[param_id] = 1u;

        const scalar theta = params[param_id].theta();
        if (theta <= 0.f || theta >= 2.f * constant<traccc::scalar>::pi) {
            TRACCC_ERROR_DEVICE("Theta is zero after propagation");
            params_liveness[param_id] = 0u;
        }

        if (!std::isfinite(params[param_id].phi())) {
            TRACCC_ERROR_DEVICE(
                "Phi is infinite after propagation (Matrix inversion)");
            params_liveness[param_id] = 0u;
        }

        if (math::fabs(params[param_id].qop()) == 0.f) {
            TRACCC_ERROR_DEVICE("q/p is zero after propagation");
            params_liveness[param_id] = 0u;
        }
    } else {
        params_liveness[param_id] = 0u;
    }

    if (params_liveness[param_id] == 0 &&
        n_cands >= cfg.min_track_candidates_per_track) {
        TRACCC_VERBOSE_DEVICE("Create tip: No next sensitive found");
        auto tip_pos = tips.push_back(link_idx);
        tip_lengths.at(tip_pos) = n_cands;
    }
}

}  // namespace traccc::device
