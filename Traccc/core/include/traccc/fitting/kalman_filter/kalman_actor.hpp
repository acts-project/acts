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
#include "traccc/finding/measurement_selector.hpp"
#include "traccc/fitting/kalman_filter/gain_matrix_updater.hpp"
#include "traccc/fitting/kalman_filter/is_line_visitor.hpp"
#include "traccc/fitting/kalman_filter/two_filters_smoother.hpp"
#include "traccc/fitting/status_codes.hpp"
#include "traccc/utils/logging.hpp"
#include "traccc/utils/particle.hpp"

// detray include(s).
#include <detray/definitions/navigation.hpp>
#include <detray/navigation/direct_navigator.hpp>
#include <detray/propagator/actors/parameter_updater.hpp>
#include <detray/propagator/actors/surface_sequencer.hpp>
#include <detray/propagator/base_actor.hpp>

// vecmem include(s)
#include <vecmem/containers/device_vector.hpp>

namespace traccc {

enum class kalman_actor_direction {
    FORWARD_ONLY,
    BACKWARD_ONLY,
    BIDIRECTIONAL
};

template <typename algebra_t, typename surface_t>
struct kalman_actor_state {

    using sequencer_t =
        typename detray::actor::surface_sequencer<surface_t>::state;

    /// Constructor with the vector of track states
    TRACCC_HOST_DEVICE
    kalman_actor_state(
        const typename edm::track_collection<algebra_t>::device::proxy_type&
            track,
        const typename edm::track_state_collection<algebra_t>::device&
            track_states,
        const edm::measurement_collection::const_device& measurements,
        vecmem::device_vector<surface_t> sequence,
        const measurement_selector::config& calib_cfg)
        : m_calib_cfg{calib_cfg},
          m_track{track},
          m_track_states{track_states},
          m_measurements{measurements},
          m_sequencer{sequence} {

        reset();
    }

    /// Get the track state at a given position along the track
    TRACCC_HOST_DEVICE
    typename edm::track_state_collection<algebra_t>::device::proxy_type at(
        unsigned int i) {
        assert(m_track.constituent_links().at(i).type ==
               edm::track_constituent_link::track_state);
        return m_track_states.at(m_track.constituent_links().at(i).index);
    }

    /// Get the track state at a given position along the track
    TRACCC_HOST_DEVICE
    typename edm::track_state_collection<algebra_t>::device::const_proxy_type
    at(unsigned int i) const {
        assert(m_track.constituent_links().at(i).type ==
               edm::track_constituent_link::track_state);
        return m_track_states.at(m_track.constituent_links().at(i).index);
    }

    /// @return the reference of track state currently handled by the KF
    TRACCC_HOST_DEVICE
    typename edm::track_state_collection<algebra_t>::device::proxy_type
    operator()(const bool backward_mode) {
        const auto state_idx{
            backward_mode ? static_cast<int>(size() - 1u) - m_idx : m_idx};
        assert(state_idx >= 0);
        assert(state_idx < static_cast<int>(size()));
        return at(static_cast<unsigned int>(state_idx));
    }

    /// Reset the state
    TRACCC_HOST_DEVICE
    void reset() {
        m_idx = 0;
        n_holes = 0u;
    }

    /// Advance the next track state
    TRACCC_HOST_DEVICE
    void next() { m_idx++; }

    /// @TODO: Const-correctness broken due to a vecmem bug
    /// @returns the number of track states
    TRACCC_HOST_DEVICE
    unsigned int size() /*const*/ { return m_track.constituent_links().size(); }

    /// @return access to the surface sequencer - non-const
    TRACCC_HOST_DEVICE sequencer_t& sequencer() { return m_sequencer; }

    /// @return access to the surface sequencer - non-const
    TRACCC_HOST_DEVICE const sequencer_t& sequencer() const {
        return m_sequencer;
    }

    /// Add the next surface to the sequence
    TRACCC_HOST_DEVICE void add_to_sequence(
        typename sequencer_t::surface_type sf_desc) {
        assert(!sf_desc.identifier().is_invalid());

        // Don't collect data
        if (m_sequencer.sequence().capacity() == 0u) {
            return;
        }

        m_sequencer.sequence().push_back(sf_desc);
        DETRAY_VERBOSE_HOST("Added: " << sf_desc);
    }

    /// @return true if the we reach the end of the track state vector
    TRACCC_HOST_DEVICE
    bool finished() /*const*/ { return m_idx == static_cast<int>(size()); }

    /// @TODO: Const-correctness broken due to a vecmem bug
    /// @returns the current number of missed states during forward fit
    TRACCC_HOST_DEVICE
    unsigned int count_missed_fit() /*const*/ {
        unsigned int n_missed{0u};

        for (unsigned int i = 0u; i < size(); ++i) {
            const edm::track_state trk_state = at(i);
            if (!trk_state.is_hole() &&
                trk_state.filtered_params().is_invalid()) {
                TRACCC_DEBUG_HOST_DEVICE(
                    "Missed track state %d/%d on surface %d during forward fit",
                    i, size(), at(i).filtered_params().surface_link().index());
                ++n_missed;
            }
        }

        return n_missed;
    }

    /// @TODO: Const-correctness broken due to a vecmem bug
    /// @returns the current number of missed states during smoothing
    TRACCC_HOST_DEVICE
    unsigned int count_missed_smoother() /*const*/ {
        unsigned int n_missed{0u};

        for (unsigned int i = 0u; i < size(); ++i) {
            const edm::track_state trk_state = at(i);
            if (!trk_state.is_hole() &&
                trk_state.smoothed_params().is_invalid()) {
                TRACCC_DEBUG_HOST_DEVICE(
                    "Missed track state %d/%d on surface %d during smoothing",
                    i, size(), at(i).smoothed_params().surface_link().index());
                ++n_missed;
            }
        }

        return n_missed;
    }

    /// @return true if the current surface could be matched to a track state
    /// @TODO: Remove once direct navigator is used in forward pass
    template <typename propagation_state_t>
    TRACCC_HOST_DEVICE bool match_surface_to_track_state(
        propagation_state_t& propagation) {

        // Surface matching logic not needed for direct navigation
        if constexpr (std::same_as<
                          typename propagation_state_t::navigator_state_type,
                          detray::direct_navigator<
                              typename propagation_state_t::detector_type>>) {
            return true;
        } else {
            const auto& navigation = propagation.navigation();
            const bool backward_mode{navigation.direction() ==
                                     detray::navigation::direction::e_backward};
            edm::track_state trk_state = (*this)(backward_mode);

            /// Check if a surface should be counted as a hole
            auto check_if_hole = [this]<typename nav_state_t>(
                                     const nav_state_t& nav) {
                if (this->do_precise_hole_count || !nav.current().is_edge()) {
                    TRACCC_VERBOSE_HOST_DEVICE("---> State flagged as hole");
                    ++(this->n_holes);
                    return true;
                }
                return false;
            };

            TRACCC_VERBOSE_HOST(
                "Navigation found: " << navigation.geometry_identifier());

            // Surface was matched, continue with KF algorithm if state not
            // a hole
            if (navigation.geometry_identifier() ==
                trk_state.filtered_params().surface_link()) {
                // Count a hole, if track finding did not find a measurement
                if (trk_state.is_hole()) {
                    TRACCC_VERBOSE_HOST_DEVICE(
                        "-> Track finding flagged this as hole");
                    check_if_hole(navigation);
                }

                TRACCC_VERBOSE_HOST_DEVICE(
                    "-> Matched current surface to next track state: %d/%d",
                    m_idx + 1, size());

                // If track finding did not find measurement on this
                // surface: skip
                return !trk_state.is_hole();
            }

            // Check if the current navigation surfaces can be found on a
            // later track state. That means the current track state was
            // skipped by the navigator: Advance the internal track state
            // index
            TRACCC_VERBOSE_HOST_DEVICE(
                "--> Check other states on track for a match");

            unsigned int n{1u};
            constexpr bool skip_this_surface{false};

            /// If a matching state was found, apply distance to internal
            /// state idx
            auto handle_skipped =
                [this, backward_mode](
                    const int i, const unsigned int n_skipped,
                    const detray::geometry::identifier geo_id) {
                    const int n_states{static_cast<int>(this->size())};
                    const int state_idx{backward_mode ? n_states - 1 - i : i};

                    assert(state_idx >= 0);
                    assert(state_idx < n_states);

                    if (this->at(static_cast<unsigned int>(state_idx))
                            .filtered_params()
                            .surface_link() == geo_id) {

                        TRACCC_VERBOSE_HOST_DEVICE(
                            "--> Matched to later state on track: "
                            "navigator skipped"
                            " %d surfaces in between",
                            n_skipped);

                        this->m_idx += static_cast<int>(n_skipped);
                        assert(this->m_idx < n_states);

                        return true;
                    }
                    return false;
                };

            // If the last track state is not matched, the current surface
            // must be a hole
            if (m_idx + 1 == static_cast<int>(size())) {
                TRACCC_DEBUG_HOST_DEVICE("--> Evaluate last state");
                check_if_hole(navigation);
                return skip_this_surface;
            }

            // If we are on the last state and the navigation surface does
            // not match, it must be an additional surface
            // -> continue navigation until matched
            assert(m_idx + 1 < static_cast<int>(size()));
            for (int i = m_idx + 1; i < static_cast<int>(size()); ++i) {
                // Found matching surface at track state 'i'
                if (handle_skipped(i, n, navigation.geometry_identifier())) {
                    fit_result =
                        backward_mode
                            ? kalman_fitter_status::ERROR_SMOOTHER_SKIPPED_STATE
                            : kalman_fitter_status::ERROR_UPDATER_SKIPPED_STATE;
                    return !skip_this_surface;
                }
                ++n;
            }
            assert(m_idx >= 0);
            assert(m_idx < static_cast<int>(size()));

            // Mismatch was not from missed state: Is a hole
            TRACCC_DEBUG_HOST_DEVICE(
                "--> Did NOT find state: might be hole...");
            if (check_if_hole(navigation)) {
                TRACCC_DEBUG_HOST("--> Expected surfaces (along the track):");
                for (unsigned int i = 0u; i < size(); ++i) {
                    TRACCC_DEBUG_HOST(
                        "   - " << at(i).filtered_params().surface_link());
                }
            } else if (navigation.current().is_edge()) {
                TRACCC_VERBOSE_HOST_DEVICE("--> Hit surface edge: Not a hole");
            }

            // After additional surface, keep navigating until match is found
            return skip_this_surface;
        }
    }

    /// Measurement calibration configuration
    measurement_selector::config m_calib_cfg{};

    /// Object describing the track fit
    typename edm::track_collection<algebra_t>::device::proxy_type m_track;
    /// All track states in the event
    typename edm::track_state_collection<algebra_t>::device m_track_states;
    /// All measurements in the event
    edm::measurement_collection::const_device m_measurements;

    /// The surface sequencer
    sequencer_t m_sequencer;

    /// Index of the current track state
    int m_idx;

    /// The number of holes (The number of sensitive surfaces which do not
    /// have a measurement for the track pattern)
    unsigned int n_holes{0u};

    /// Result of the fitter pass
    kalman_fitter_status fit_result = kalman_fitter_status::SUCCESS;

    /// Finish the navigation beyond the track states in the fitter to find
    /// all holes
    bool do_precise_hole_count = false;
};

/// Detray actor for Kalman filtering
template <typename algebra_t, typename surface_t,
          kalman_actor_direction direction_e>
struct kalman_actor : detray::base_actor {

    // Actor state
    using state = kalman_actor_state<algebra_t, surface_t>;

    /// Actor operation to perform the Kalman filtering
    ///
    /// @param actor_state the actor state
    /// @param propagation the propagator state
    template <typename propagator_state_t>
    TRACCC_HOST_DEVICE void operator()(
        state& actor_state, propagator_state_t& propagation,
        detray::actor::parameter_transporter_result<algebra_t>& res) const {

        auto& stepping = std::as_const(propagation).stepping();
        auto& navigation = std::as_const(propagation).navigation();

        TRACCC_VERBOSE_HOST_DEVICE("Actor: Kalman Fitter (status %d)...",
                                   actor_state.fit_result);

        // Allow to count holes after the intial track states
        if (actor_state.do_precise_hole_count && actor_state.finished()) {
            if (navigation.is_on_sensitive()) {
                TRACCC_VERBOSE_HOST_DEVICE(
                    "Track is already complete: This surface is a hole");
                // At this point every surface is a hole
                actor_state.n_holes++;
            }
            return;
        }

        // triggered only for sensitive surfaces
        if (navigation.is_on_sensitive()) {

            // Dynamically check the navigation direction
            const bool backward_mode{navigation.direction() ==
                                     detray::navigation::direction::e_backward};

            TRACCC_VERBOSE_HOST(
                "Expected: "
                << actor_state(backward_mode).filtered_params().surface_link());

            // Increase the hole count if the propagator stops at an additional
            // surface and wait for the next sensitive surface to match
            if (!actor_state.match_surface_to_track_state(propagation)) {
                // Add this to the surface sequence for the backward fit
                if (!backward_mode &&
                    navigation.current_surface().has_material()) {
                    actor_state.add_to_sequence(navigation.current().surface());
                }
                return;
            } else if (actor_state.fit_result !=
                           kalman_fitter_status::SUCCESS &&
                       (!actor_state.do_precise_hole_count ||
                        res.destination_params().is_invalid())) {

                // Make sure this is still counted as a reached surface even
                // though the fitter bailed out (as indicated by fitter status)
                if (actor_state.do_precise_hole_count &&
                    res.destination_params().is_invalid()) {
                    actor_state(backward_mode)
                        .filtered_params()
                        .set_parameter_vector(
                            navigation.current_surface().free_to_bound_vector(
                                {}, stepping()));
                    actor_state(backward_mode)
                        .filtered_params()
                        .set_covariance(matrix::identity<
                                        traccc::bound_matrix<algebra_t>>());
                    actor_state(backward_mode)
                        .filtered_params()
                        .set_surface_link(navigation.geometry_identifier());

                    actor_state(backward_mode)
                        .smoothed_params()
                        .set_parameter_vector(
                            navigation.current_surface().free_to_bound_vector(
                                {}, stepping()));
                    actor_state(backward_mode)
                        .smoothed_params()
                        .set_covariance(matrix::identity<
                                        traccc::bound_matrix<algebra_t>>());
                    actor_state(backward_mode)
                        .smoothed_params()
                        .set_surface_link(navigation.geometry_identifier());
                }
                // Surface matched but encountered error: Abort fit
                propagation.navigation().abort(
                    fitter_debug_msg{actor_state.fit_result});
                propagation.heartbeat(false);
                return;
            }

            auto& sequencer = actor_state.sequencer();
            if (sequencer.sequence().capacity() > 0u &&
                sequencer.sequence().size() ==
                    sequencer.sequence().capacity()) {
                DETRAY_ERROR_HOST_DEVICE("Sequence overflow!");
                sequencer.set_overflow();
                propagation.navigation().exit();
                propagation.heartbeat(false);
                return;
            }

            // Fetch matched track state
            edm::track_state trk_state = actor_state(backward_mode);
            bound_track_parameters<algebra_t>& bound_param =
                res.destination_params();

            assert(trk_state.filtered_params().surface_link() ==
                   bound_param.surface_link());

            // Run Kalman Gain Updater
            const auto sf = navigation.current_surface();
            const bool is_line = detail::is_line(sf);

            const auto measurement =
                actor_state.m_measurements.at(trk_state.measurement_index());

            if (!backward_mode) {
                if constexpr (direction_e ==
                                  kalman_actor_direction::FORWARD_ONLY ||
                              direction_e ==
                                  kalman_actor_direction::BIDIRECTIONAL) {
                    // Wrap the phi and theta angles in their valid ranges
                    normalize_angles<algebra_t>(bound_param);

                    // Forward filter
                    TRACCC_DEBUG_HOST_DEVICE("Run filtering...");
                    actor_state.fit_result = gain_matrix_updater<algebra_t>{}(
                        trk_state, measurement, bound_param,
                        actor_state.m_calib_cfg, is_line);

                    // Update the propagation flow
                    bound_param = trk_state.filtered_params();

                    // Add this to the surface sequence for the backward fit
                    actor_state.add_to_sequence(navigation.current().surface());

                    if (actor_state.fit_result ==
                        kalman_fitter_status::SUCCESS) {

                        // Calculate the chi2 on the filtered parameters
                        trk_state.filtered_chi2() =
                            measurement_selector::predicted_chi2(
                                measurement, bound_param,
                                actor_state.m_calib_cfg, is_line);
                    }
                } else {
                    assert(false);
                }
            } else {
                if constexpr (direction_e ==
                                  kalman_actor_direction::BACKWARD_ONLY ||
                              direction_e ==
                                  kalman_actor_direction::BIDIRECTIONAL) {
                    // Backward filter for smoothing
                    TRACCC_DEBUG_HOST_DEVICE("Run smoothing...");

                    // Forward filter did not find this state: cannot
                    // smoothe
                    if (trk_state.filtered_params().is_invalid()) {
                        TRACCC_ERROR_HOST_DEVICE(
                            "Track state not filtered by forward fit. "
                            "Skipping");
                        actor_state.fit_result =
                            kalman_fitter_status::ERROR_UPDATER_SKIPPED_STATE;
                    } else {
                        actor_state.fit_result =
                            two_filters_smoother<algebra_t>{}(
                                trk_state, measurement, bound_param,
                                actor_state.m_calib_cfg, is_line);
                    }
                } else {
                    assert(false);
                }
            }

            // Abort if the Kalman update fails
            if (actor_state.fit_result != kalman_fitter_status::SUCCESS) {
                TRACCC_ERROR_HOST_DEVICE(
                    "Fit failure: KF status %d",
                    static_cast<int>(actor_state.fit_result));
                TRACCC_ERROR_HOST(" --> reason: " << fitter_debug_msg{
                                      actor_state.fit_result}());

                if (actor_state.do_precise_hole_count ||
                    bound_param.is_invalid()) {
                    if (bound_param.is_invalid()) {
                        trk_state.filtered_params().set_parameter_vector(
                            navigation.current_surface().free_to_bound_vector(
                                {}, stepping()));
                        trk_state.filtered_params().set_covariance(
                            matrix::identity<
                                traccc::bound_matrix<algebra_t>>());
                        trk_state.filtered_params().set_surface_link(
                            navigation.geometry_identifier());

                        trk_state.smoothed_params().set_parameter_vector(
                            navigation.current_surface().free_to_bound_vector(
                                {}, stepping()));
                        trk_state.smoothed_params().set_covariance(
                            matrix::identity<
                                traccc::bound_matrix<algebra_t>>());
                        trk_state.smoothed_params().set_surface_link(
                            navigation.geometry_identifier());
                    }
                    propagation.navigation().abort(
                        fitter_debug_msg{actor_state.fit_result});
                    propagation.heartbeat(false);
                    return;
                }
            }

            // Change the charge of hypothesized particles when the sign of qop
            // is changed (This rarely happens when qop is set with a poor seed
            // resolution)
            propagation.set_particle(detail::correct_particle_hypothesis(
                stepping.particle_hypothesis(), bound_param));

            // Update the expected track state
            actor_state.next();

            // No need to continue (unless need precise hole count)
            if (actor_state.finished() && !actor_state.do_precise_hole_count) {
                TRACCC_VERBOSE_HOST_DEVICE("Kalman Actor: finished");
                propagation.navigation().exit();
                propagation.heartbeat(false);
                return;
            }

            // Flag renavigation of the current candidate (unless for
            // overlap)
            if (math::fabs(navigation()) > 1.f * unit<float>::um) {
                propagation.navigation().set_high_trust();
            } else {
                TRACCC_DEBUG_HOST_DEVICE(
                    "Encountered overlap, jump to next surface");
            }

            // Signal that paramter update is needed
            res.status = detray::actor::status::e_success;
        } else if (navigation.direction() ==
                       detray::navigation::direction::e_forward &&
                   navigation.encountered_sf_material()) {
            assert(navigation.is_on_passive() || navigation.is_on_portal());

            // Add this to the surface sequence for the backward fit
            actor_state.add_to_sequence(navigation.current().surface());
        }
    }
};

}  // namespace traccc
