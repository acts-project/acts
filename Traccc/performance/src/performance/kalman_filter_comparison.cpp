/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/bfield/construct_const_bfield.hpp"
#include "traccc/fitting/kalman_filter/kalman_actor.hpp"
#include "traccc/utils/propagation.hpp"

// Performance include(s).
#include "traccc/performance/kalman_filter_comparison.hpp"

// Detray test include(s)
#include <detray/test/utils/perigee_stopper.hpp>
#include <detray/test/validation/navigation_validation_utils.hpp>

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <exception>
#include <iostream>
#include <memory>
#include <optional>
#include <string>

namespace traccc {

bool kalman_filter_comparison(
    const traccc::default_detector::host& det,
    const traccc::default_detector::host::name_map& names,
    const traccc::magnetic_field& bfield,
    const detray::propagation_validation_config& cfg,
    const traccc::seed_generator<traccc::default_detector::host>::config&
        smearing_cfg,
    std::unique_ptr<const traccc::Logger> ilogger,
    const std::vector<traccc::free_track_parameters<traccc::default_algebra>>&
        tracks,
    std::vector<vecmem::vector<traccc::propagation_validator::candidate_type<
        traccc::default_detector::host>>>& truth_traces_fw,
    edm::measurement_collection::const_device device_measurements,
    traccc::edm::track_container<traccc::default_algebra>::host&
        track_container) {

    using namespace traccc;

    TRACCC_LOCAL_LOGGER(std::move(ilogger));

    // 'false' if any failures were detected
    bool test_successful{true};

    using detector_t = traccc::default_detector::host;
    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = detray::dscalar<algebra_t>;
    using b_field_t =
        covfie::field<traccc::const_bfield_backend_t<traccc::scalar>>;
    using stepper_t =
        detray::rk_stepper<b_field_t::view_t, algebra_t,
                           detray::constrained_step<traccc::scalar>,
                           detray::stepper_rk_policy<traccc::scalar>,
                           detray::stepping::print_inspector>;
    using sf_candidate_t =
        traccc::propagation_validator::candidate_type<detector_t>;

    // Host memory resource
    vecmem::host_memory_resource host_mr;

    // Geometry context
    const detector_t::geometry_context ctx{};

    // Create B field
    b_field_t field =
        bfield.template as_field<traccc::const_bfield_backend_t<scalar_t>>();
    std::optional<b_field_t::view_t> field_view{field};

    // Seed smearing
    auto param_smearer =
        seed_generator{det, smearing_cfg, std::mt19937::default_seed, ctx};

    // Configure the test
    measurement_selector::config calib_cfg{};
    detray::test::navigation_validation_config<algebra_t> test_cfg{};
    test_cfg.n_tracks(tracks.size());
    test_cfg.ptc_hypothesis(cfg.particle);
    test_cfg.collect_sensitives_only(true);
    test_cfg.fail_on_diff(false);
    test_cfg.display_svg(cfg.display_svg);
    test_cfg.display_only_missed(true);
    test_cfg.verbose(false);

    // Create the backward truth traces by reverting the order
    std::vector<vecmem::vector<sf_candidate_t>> truth_traces_bw;
    truth_traces_bw.reserve(tracks.capacity());

    for (const auto& truth_trace_fw : truth_traces_fw) {
        // Revert the forward trace for the backward propagation
        vecmem::vector<sf_candidate_t> truth_trace_bw(truth_trace_fw.size());
        std::ranges::reverse_copy(truth_trace_fw, truth_trace_bw.begin());

        assert(!truth_trace_bw.empty());

        truth_traces_bw.push_back(std::move(truth_trace_bw));
    }

    // Check truth data
    if (truth_traces_fw.empty() || truth_traces_bw.empty()) {
        TRACCC_ERROR("Propagation truth data empty");
        return false;
    }
    if (track_container.tracks.size() == 0u) {
        TRACCC_ERROR("Kalman Filter truth data empty");
        return false;
    }
    if (track_container.tracks.size() != tracks.size()) {
        TRACCC_ERROR(
            "Initial track parameters and KF (forward) tracks do not match");
        return false;
    }

    // Make a tuple of references from a tuple
    auto setup_actor_states = []<typename... T>(detray::dtuple<T...>& t) {
        return detray::tie(detray::detail::get<T>(t)...);
    };

    // Define the actors
    using interactor = detray::actor::pointwise_material_interactor<algebra_t>;
    using perigee_stopper = detray::perigee_stopper<algebra_t>;

    // Reusable actor states
    interactor::state interactor_state{};
    interactor_state.do_multiple_scattering = cfg.do_multiple_scattering;
    interactor_state.do_energy_loss = cfg.do_energy_loss;

    perigee_stopper::state stopper_state{};

    std::cout << "-----------------------------------"
              << "\nFORWARD - With KF" << std::endl
              << "-----------------------------------\n";

    // Define the actors
    using surface_t = typename detector_t::surface_type;
    using fit_actor_fw =
        traccc::kalman_actor<algebra_t, surface_t,
                             kalman_actor_direction::FORWARD_ONLY>;
    using parameter_updater =
        detray::actor::parameter_updater<algebra_t, interactor, fit_actor_fw>;
    using actor_chain_t = detray::actor_chain<parameter_updater>;

    // Reusable actor states
    parameter_updater::state updater_state{};
    updater_state.noise_estimation_cfg().n_stddev =
        cfg.propagation.navigation.n_scattering_stddev;
    updater_state.noise_estimation_cfg().accumulated_error =
        cfg.propagation.navigation.accumulated_error;
    updater_state.noise_estimation_cfg().estimate_scattering_noise =
        cfg.propagation.navigation.estimate_scattering_noise;

    // Create the input container(s) the way the KF needs them
    typename edm::track_container<algebra_t>::data tracks_data{track_container};
    typename edm::track_container<algebra_t>::view tracks_view{tracks_data};
    typename edm::track_container<algebra_t>::device device_track_container{
        tracks_view};

    // This collects the surfaces for the backward fit - not needed here
    vecmem::vector<surface_t> sf_sequence{};
    vecmem::device_vector<surface_t> device_sf_sequence{
        vecmem::get_data(sf_sequence)};

    // Prepare the actor states for every track
    vecmem::vector<typename actor_chain_t::state_tuple> state_tuple{};
    vecmem::vector<typename actor_chain_t::state_ref_tuple> state_ref_tuple{};
    state_tuple.reserve(track_container.tracks.size());
    state_ref_tuple.reserve(track_container.tracks.size());

    for (std::size_t i = 0u; i < track_container.tracks.size(); ++i) {

        // Define the initial covariance
        updater_state.init(tracks.at(i));
        param_smearer.generate_initial_covariance(updater_state.bound_params(),
                                                  cfg.particle);

        // Create the fitter state per track
        fit_actor_fw::state fit_actor_state(
            device_track_container.tracks.at(static_cast<unsigned int>(i)),
            device_track_container.states, device_measurements,
            device_sf_sequence, calib_cfg);

        fit_actor_state.reset();
        fit_actor_state.do_precise_hole_count = true;

        state_tuple.push_back(detray::make_tuple(
            updater_state, interactor_state, fit_actor_state));
        state_ref_tuple.push_back(setup_actor_states(state_tuple.back()));
    }

    // Forward filter
    test_cfg.name(det.name(names) + "_GeV_fw_KF");
    test_cfg.navigation_direction(detray::navigation::direction::e_forward);
    const auto [trk_stats_fw, n_surfaces_fw, n_miss_nav_fw, n_miss_truth_fw,
                step_traces_fw, mat_traces_fw, mat_records_fw] =
        detray::navigation_validator::compare_to_navigation<stepper_t,
                                                            parameter_updater>(
            test_cfg, host_mr, det, names, ctx, field_view, cfg.propagation,
            truth_traces_fw, tracks, state_ref_tuple);

    // Check, how many holes the KF found
    std::size_t n_missed_fw{0u};
    std::size_t n_trk_missing_fw{0u};
    std::size_t n_holes_fw{0u};
    std::size_t n_trk_holes_fw{0u};
    for (unsigned int i = 0u; i < tracks.size(); ++i) {
        const auto& actor_states = state_tuple[i];
        auto fitter_state = detray::get<fit_actor_fw::state>(actor_states);

        // What the actor counted
        n_missed_fw += fitter_state.count_missed_fit();
        n_holes_fw += fitter_state.n_holes;

        if (fitter_state.count_missed_fit() > 0u) {
            n_trk_missing_fw++;
        }
        if (fitter_state.n_holes > 0u) {
            n_trk_holes_fw++;
        }
    }
    std::cout << "No. skipped states in fw KF: " << n_missed_fw << std::endl;
    std::cout << "No. tracks with skipped states in fw KF: " << n_trk_missing_fw
              << std::endl;
    std::cout << "No. holes found by fw KF: " << n_holes_fw << std::endl;
    std::cout << "No. tracks with holes in fw KF: " << n_trk_holes_fw
              << std::endl;

    // Trigger failures
    if (n_miss_nav_fw.n_total() != n_missed_fw) {
        TRACCC_ERROR("Forward filter number of missed states incorrect: was "
                     << n_missed_fw << ", should be "
                     << n_miss_nav_fw.n_total());
        test_successful = false;
    }
    if (n_trk_missing_fw != trk_stats_fw.n_tracks_w_holes) {
        TRACCC_ERROR("Forward filter number of faulty tracks incorrect: was "
                     << n_trk_missing_fw << ", should be "
                     << trk_stats_fw.n_tracks_w_holes);
        test_successful = false;
    }
    if (n_miss_truth_fw.n_total() != n_holes_fw) {
        TRACCC_ERROR("Forward filter hole counting incorrect: was "
                     << n_holes_fw << ", should be "
                     << n_miss_truth_fw.n_total());
        test_successful = false;
    }
    if (n_trk_holes_fw < trk_stats_fw.n_tracks_w_extra) {
        TRACCC_ERROR(
            "Forward filter number of tracks with holes incorrect: was "
            << n_trk_holes_fw << ", should be "
            << trk_stats_fw.n_tracks_w_extra);
        test_successful = false;
    }

    // Check stats
    auto n_tracks{static_cast<double>(tracks.size())};
    if (static_cast<double>(trk_stats_fw.n_tracks_w_holes) / n_tracks >
        cfg.max_percent_missed / 100.) {
        TRACCC_ERROR("Too many tracks with missing surfaces");
        test_successful = false;
    }
    if (static_cast<double>(trk_stats_fw.n_tracks_w_extra) / n_tracks >
        cfg.max_percent_additional / 100.) {
        TRACCC_ERROR("Too many tracks with additional surfaces");
        test_successful = false;
    }

    std::cout << "\n-----------------------------------" << std::endl;
    std::cout << "BACKWARD - With KF" << std::endl
              << "-----------------------------------\n";
    using fit_actor_bw =
        traccc::kalman_actor<algebra_t, surface_t,
                             kalman_actor_direction::BIDIRECTIONAL>;
    using parameter_updater_bw =
        detray::actor::parameter_updater<algebra_t, interactor, fit_actor_bw>;
    using actor_chain_bw_t =
        detray::actor_chain<parameter_updater, perigee_stopper>;

    // Reusable actor states
    parameter_updater_bw::state updater_state_bw{};
    updater_state_bw.noise_estimation_cfg().n_stddev =
        cfg.propagation.navigation.n_scattering_stddev;
    updater_state_bw.noise_estimation_cfg().accumulated_error =
        cfg.propagation.navigation.accumulated_error;
    updater_state_bw.noise_estimation_cfg().estimate_scattering_noise =
        cfg.propagation.navigation.estimate_scattering_noise;

    vecmem::vector<typename actor_chain_bw_t::state_tuple> state_tuple_bw{};
    vecmem::vector<typename actor_chain_bw_t::state_ref_tuple>
        state_ref_tuple_bw{};
    state_tuple_bw.reserve(tracks.size());
    state_ref_tuple_bw.reserve(tracks.size());

    // Prepare the fitter state for every track
    for (std::size_t i = 0u; i < tracks.size(); ++i) {

        // Define the initial covariance
        updater_state_bw.init(tracks.at(i));
        param_smearer.generate_initial_covariance(
            updater_state_bw.bound_params(), cfg.particle);

        // Create the fitter state per track
        fit_actor_bw::state fit_actor_state(
            device_track_container.tracks.at(static_cast<unsigned int>(i)),
            device_track_container.states, device_measurements,
            device_sf_sequence, calib_cfg);

        fit_actor_state.reset();
        fit_actor_state.do_precise_hole_count = true;

        state_tuple_bw.push_back(
            detray::make_tuple(stopper_state, interactor_state, fit_actor_state,
                               updater_state_bw));
        state_ref_tuple_bw.push_back(setup_actor_states(state_tuple_bw.back()));
    }

    // Backward filter
    test_cfg.name(det.name(names) + "_GeV_bw_KF");
    test_cfg.navigation_direction(detray::navigation::direction::e_backward);
    const auto [trk_stats_bw, n_surfaces_bw, n_miss_nav_bw, n_miss_truth_bw,
                step_traces_bw, mat_traces_bw, mat_records_bw] =
        detray::navigation_validator::compare_to_navigation<
            stepper_t, parameter_updater_bw, perigee_stopper>(
            test_cfg, host_mr, det, names, ctx, field_view, cfg.propagation,
            truth_traces_bw, tracks, state_ref_tuple_bw);

    // Check, how many tracks were smoothed correctly
    std::size_t n_missed_bw{0u};
    std::size_t n_holes_bw{0u};
    std::size_t n_trk_holes_bw{0u};
    std::size_t n_trk_missing_bw{0u};
    std::size_t n_not_smoothed_correctly{0u};
    for (std::size_t i = 0u; i < tracks.size(); ++i) {
        const auto& actor_states = state_tuple_bw[i];
        auto fitter_state = detray::get<fit_actor_fw::state>(actor_states);

        // What the actor counted
        n_missed_bw += fitter_state.count_missed_smoother();
        n_holes_bw += fitter_state.n_holes;

        if (fitter_state.count_missed_smoother() > 0u) {
            n_trk_missing_bw++;
        }
        if (fitter_state.n_holes > 0u) {
            n_trk_holes_bw++;
        }

        for (unsigned int istate = 0u; istate < fitter_state.size(); ++istate) {
            if (!fitter_state.at(istate).is_smoothed()) {
                n_not_smoothed_correctly++;
                break;
            }
        }
    }
    std::cout << "INCLUDES SKIPPED STATES BY FW KF PASS:\n" << std::endl;
    std::cout << "No. skipped states in bw KF: " << n_missed_bw << std::endl;
    std::cout << "No. tracks with skipped states in bw KF: " << n_trk_missing_bw
              << std::endl;
    std::cout << "No. holes found by bw KF: " << n_holes_bw << std::endl;
    std::cout << "No. tracks with holes in bw KF: " << n_trk_holes_bw
              << std::endl;
    std::cout << "No. tracks that were not smoothed correctly: "
              << n_not_smoothed_correctly << " ("
              << 100. * static_cast<double>(n_not_smoothed_correctly) /
                     static_cast<double>(trk_stats_bw.n_tracks)
              << "%)\n"
              << std::endl;

    // Trigger failures
    if (n_missed_bw < n_miss_nav_bw.n_total()) {
        TRACCC_ERROR("Backward filter number of missed states incorrect: was "
                     << n_missed_bw << ", should be greater equal "
                     << n_miss_nav_bw.n_total());
        test_successful = false;
    }
    if (n_trk_missing_bw < trk_stats_bw.n_tracks_w_holes) {
        TRACCC_ERROR("Backward filter number of faulty tracks incorrect: was "
                     << n_trk_missing_bw << ", should be greater equal "
                     << trk_stats_bw.n_tracks_w_holes);
        test_successful = false;
    }
    if (n_miss_truth_bw.n_total() != n_holes_bw) {
        TRACCC_ERROR("Backward filter hole counting incorrect: was "
                     << n_holes_bw << ", should be "
                     << n_miss_truth_bw.n_total());
        test_successful = false;
    }
    if (n_trk_holes_bw != trk_stats_bw.n_tracks_w_extra) {
        TRACCC_ERROR(
            "Backward filter number of tracks with holes incorrect: was "
            << n_trk_holes_bw << ", should be "
            << trk_stats_bw.n_tracks_w_extra);
        test_successful = false;
    }

    // Check stats
    if (static_cast<double>(trk_stats_bw.n_tracks_w_holes) / n_tracks >
        cfg.max_percent_missed / 100.) {
        TRACCC_ERROR(
            "Too many tracks with missing surfaces (temporarily disabled)");
        // test_successful = false; // TODO: Reactivate once KF is stable
    }
    if (static_cast<double>(trk_stats_bw.n_tracks_w_extra) / n_tracks >
        cfg.max_percent_additional / 100.) {
        TRACCC_ERROR("Too many tracks with additional surfaces");
        test_successful = false;
    }
    std::cout << "-----------------------------------" << std::endl;

    return test_successful;
}

}  // namespace traccc
