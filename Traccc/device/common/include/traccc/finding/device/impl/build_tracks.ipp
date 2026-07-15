/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include <type_traits>

#include "traccc/edm/measurement_helpers.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/edm/track_state_helpers.hpp"
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/utils/matrix_helpers.hpp"
#include "traccc/utils/prob.hpp"
#include "traccc/utils/subspace.hpp"

namespace traccc::device {

TRACCC_HOST_DEVICE inline void build_tracks(
    const global_index_t globalIndex, bool run_mbf,
    const measurement_selector::config calib_cfg,
    const build_tracks_payload& payload) {

    const edm::measurement_collection::const_device measurements(
        payload.tracks_view.measurements);

    const bound_track_parameters_collection_types::const_device seeds(
        payload.seeds_view);

    const vecmem::device_vector<const candidate_link> links(payload.links_view);

    const vecmem::device_vector<const unsigned int> tips(payload.tips_view);

    edm::track_collection<default_algebra>::device track_candidates(
        payload.tracks_view.tracks);
    edm::track_state_collection<default_algebra>::device track_states(
        payload.tracks_view.states);
    bound_track_parameters_collection_types::device link_predicted_params(
        payload.link_predicted_parameter_view);
    bound_track_parameters_collection_types::device link_filtered_params(
        payload.link_filtered_parameter_view);

    if (globalIndex >= tips.size()) {
        return;
    }

    const vecmem::device_vector<const unsigned int> tip_to_output_map(
        payload.tip_to_output_map);
    const unsigned int output_idx = tip_to_output_map.size()
                                        ? tip_to_output_map.at(globalIndex)
                                        : globalIndex;

    if (output_idx == std::numeric_limits<unsigned int>::max()) {
        return;
    }

    const unsigned int tip = tips.at(globalIndex);
    edm::track track = track_candidates.at(output_idx);

    // Get the link corresponding to tip
    unsigned int link_idx = tip;
    candidate_link L = links.at(link_idx);

    const unsigned int n_meas = measurements.size();

    // Track summary variables
    scalar ndf_sum = 0.f;
    scalar chi2_sum = 0.f;

    bound_matrix<default_algebra> big_lambda_tilde =
        matrix::zero<bound_matrix<default_algebra>>();
    bound_vector<default_algebra> small_lambda_tilde =
        matrix::zero<bound_vector<default_algebra>>();
    bound_matrix<default_algebra> accumulated_jacobian =
        matrix::identity<bound_matrix<default_algebra>>();

    bool mbf_first_processed = false;

    // Reversely iterate to fill the track candidates
    for (auto it = track.constituent_links().rbegin();
         it != track.constituent_links().rend(); it++) {

        while (L.meas_idx >= n_meas && L.step != 0u) {
            if (run_mbf) {
                accumulated_jacobian =
                    accumulated_jacobian * payload.jacobian_ptr[link_idx];
            }

            link_idx = L.previous_candidate_idx;
            L = links.at(link_idx);
        }

        assert(L.meas_idx < n_meas);

        if (run_mbf) {
            const unsigned int track_state_index =
                track_states.push_back(edm::make_track_state<default_algebra>(
                    measurements, L.meas_idx));
            auto track_state = track_states.at(track_state_index);

            track_state.set_hole(false);
            track_state.set_smoothed(true);

            // TODO: The fact that we store the chi2 three times is nonsense.
            track_state.filtered_chi2() = L.chi2;
            track_state.smoothed_chi2() = L.chi2;
            track_state.backward_chi2() = L.chi2;

            const bound_track_parameters<>& predicted_params =
                link_predicted_params.at(link_idx);
            const bound_track_parameters<>& filtered_params =
                link_filtered_params.at(link_idx);

            const bound_track_parameters<>::vector_type& predicted_vec =
                predicted_params.vector();
            const bound_track_parameters<>::covariance_type&
                predicted_covariance = predicted_params.covariance();

            bound_vector<default_algebra> small_lambda_hat;
            bound_matrix<default_algebra> big_lambda_hat;

            if (!mbf_first_processed) {
                small_lambda_hat =
                    matrix::zero<bound_vector<default_algebra>>();
                big_lambda_hat = matrix::zero<bound_matrix<default_algebra>>();
                mbf_first_processed = true;
            } else {
                small_lambda_hat = matrix::transpose(accumulated_jacobian) *
                                   small_lambda_tilde;
                big_lambda_hat = matrix::transpose(accumulated_jacobian) *
                                 big_lambda_tilde * accumulated_jacobian;
            }

            const edm::measurement meas = measurements.at(L.meas_idx);

            // Measurement data on surface
            const detray::dmatrix<default_algebra, 2, 1> meas_local =
                measurement_selector::calibrated_measurement_position<
                    default_algebra, 2>(meas, calib_cfg);

            // Spatial resolution (Measurement covariance)
            const detray::dmatrix<default_algebra, 2, 2> V =
                measurement_selector::calibrated_measurement_covariance<
                    default_algebra, 2>(meas, calib_cfg);

            // TODO: Does not work for line surfaces!
            const detray::dmatrix<default_algebra, 2, e_bound_size> H =
                measurement_selector::observation_model<default_algebra, 2>(
                    meas, predicted_params, false);

            const auto S_inv = masked_inverse<default_algebra>(
                H * predicted_covariance * matrix::transpose(H) + V,
                meas.dimensions());

            const auto K = predicted_covariance * matrix::transpose(H) * S_inv;
            const auto C_hat =
                matrix::identity<detray::dmatrix<default_algebra, e_bound_size,
                                                 e_bound_size>>() -
                K * H;
            const auto y = meas_local - H * predicted_vec;

            small_lambda_tilde = -1.f * matrix::transpose(H) * S_inv * y +
                                 matrix::transpose(C_hat) * small_lambda_hat;
            big_lambda_tilde =
                matrix::transpose(H) * S_inv * H +
                matrix::transpose(C_hat) * big_lambda_hat * C_hat;

            track_state.smoothed_params().vector() =
                predicted_vec - predicted_covariance * small_lambda_tilde;
            track_state.smoothed_params().covariance() =
                predicted_covariance -
                (predicted_covariance * big_lambda_tilde *
                 predicted_covariance);

            track_state.filtered_params() = filtered_params;

            accumulated_jacobian = payload.jacobian_ptr[link_idx];

            *it = {edm::track_constituent_link::track_state, track_state_index};
        } else {
            *it = {edm::track_constituent_link::measurement, L.meas_idx};
        }

        // Sanity check on chi2
        assert(L.chi2 < std::numeric_limits<traccc::scalar>::max());
        assert(L.chi2 >= 0.f);

        ndf_sum +=
            static_cast<scalar>(measurements.at(L.meas_idx).dimensions());
        chi2_sum += L.chi2;

        // Break the loop if the iterator is at the first candidate and fill the
        // seed and track quality
        if (it == track.constituent_links().rend() - 1) {
            if (run_mbf) {
                track.fit_outcome() = track_fit_outcome::SUCCESS;
                track.params() = track_states.at(it->index).smoothed_params();
            } else {
                track.fit_outcome() = track_fit_outcome::UNKNOWN;
                track.params() = seeds.at(L.seed_idx);
            }
            track.ndf() = ndf_sum - 5.f;
            track.chi2() = chi2_sum;
            track.pval() = prob(track.chi2(), track.ndf());
            track.nholes() = L.n_skipped;
        } else {
            link_idx = L.previous_candidate_idx;
            L = links.at(link_idx);
        }
    }

#ifndef NDEBUG
    // Assert that we did not make any duplicate track states.
    for (const auto& i : track.constituent_links()) {
        for (const auto& j : track.constituent_links()) {
            if (i.index != j.index) {
                // TODO: Re-enable me!
                // assert(measurements.at(i->index).identifier() !=
                //       measurement.at(j->index).identifier());
            }
        }
    }
#endif
}

}  // namespace traccc::device
