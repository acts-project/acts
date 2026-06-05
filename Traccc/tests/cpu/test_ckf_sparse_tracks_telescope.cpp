/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/bfield/construct_const_bfield.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/resolution/fitting_performance_writer.hpp"
#include "traccc/simulation/event_generators.hpp"
#include "traccc/simulation/simulator.hpp"
#include "traccc/utils/ranges.hpp"

// Test include(s).
#include "tests/ckf_telescope_test.hpp"
#include "traccc/utils/seed_generator.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <filesystem>
#include <string>

using namespace traccc;
// This defines the local frame test suite
TEST_P(CkfSparseTrackTelescopeTests, Run) {

    // Get the parameters
    const std::string name = std::get<0>(GetParam());
    const std::array<scalar, 3u> origin = std::get<1>(GetParam());
    const std::array<scalar, 3u> origin_stddev = std::get<2>(GetParam());
    const std::array<scalar, 2u> mom_range = std::get<3>(GetParam());
    const std::array<scalar, 2u> eta_range = std::get<4>(GetParam());
    const std::array<scalar, 2u> theta_range = eta_to_theta_range(eta_range);
    const std::array<scalar, 2u> phi_range = std::get<5>(GetParam());
    const traccc::pdg_particle<scalar> ptc = std::get<6>(GetParam());
    const unsigned int n_truth_tracks = std::get<7>(GetParam());
    const unsigned int n_events = std::get<8>(GetParam());
    const bool random_charge = std::get<9>(GetParam());

    // Performance writer
    traccc::fitting_performance_writer::config fit_writer_cfg;
    fit_writer_cfg.file_path = "performance_track_fitting_" + name + ".root";
    traccc::fitting_performance_writer fit_performance_writer(
        fit_writer_cfg, traccc::getDefaultLogger("FittingPerformanceWriter",
                                                 traccc::Logging::Level::INFO));

    /*****************************
     * Build a telescope geometry
     *****************************/

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;
    // Copy obejct
    vecmem::copy copy;

    // Read back detector file
    const std::string path = name + "/";
    traccc::host_detector detector;
    traccc::io::read_detector(
        detector, host_mr,
        std::filesystem::absolute(
            std::filesystem::path(path + "telescope_detector_geometry.json"))
            .native(),
        std::filesystem::absolute(
            std::filesystem::path(
                path + "telescope_detector_homogeneous_material.json"))
            .native());

    auto field = traccc::construct_const_bfield(std::get<13>(GetParam()));

    /***************************
     * Generate simulation data
     ***************************/

    // Track generator
    using generator_type =
        detray::random_track_generator<traccc::free_track_parameters<>,
                                       uniform_gen_t>;
    generator_type::configuration gen_cfg{};
    gen_cfg.n_tracks(n_truth_tracks);
    gen_cfg.origin(origin);
    gen_cfg.origin_stddev(origin_stddev);
    gen_cfg.phi_range(phi_range[0], phi_range[1]);
    gen_cfg.theta_range(theta_range[0], theta_range[1]);
    gen_cfg.mom_range(mom_range[0], mom_range[1]);
    gen_cfg.randomize_charge(random_charge);
    generator_type generator(gen_cfg);

    // Smearing value for measurements
    traccc::measurement_smearer<traccc::default_algebra> meas_smearer(
        smearing[0], smearing[1]);

    using writer_type = traccc::smearing_writer<
        traccc::measurement_smearer<traccc::default_algebra>>;

    typename writer_type::config smearer_writer_cfg{meas_smearer};

    // Run simulator
    const std::string full_path = io::data_directory() + path;
    std::filesystem::create_directories(full_path);
    auto sim = traccc::simulator<host_detector_type, b_field_t, generator_type,
                                 writer_type>(
        ptc, n_events, detector.as<detector_traits>(),
        field.as_field<traccc::const_bfield_backend_t<traccc::scalar>>(),
        std::move(generator), std::move(smearer_writer_cfg), full_path);

    sim.run();

    /*****************************
     * Do the reconstruction
     *****************************/

    // Seed generator
    seed_generator<host_detector_type> sg(detector.as<detector_traits>(),
                                          stddevs);

    // Finding algorithm configuration
    typename traccc::finding_config cfg;
    cfg.ptc_hypothesis = ptc;
    cfg.chi2_max = 200.f;
    cfg.min_p = 0;
    cfg.min_pT = 10.f * unit<float>::MeV;

    // Finding algorithm object
    traccc::host::combinatorial_kalman_filter_algorithm host_finding(cfg,
                                                                     host_mr);

    // Fitting algorithm object
    traccc::fitting_config fit_cfg;
    fit_cfg.ptc_hypothesis = ptc;
    fit_cfg.min_pT = 100.f * traccc::unit<float>::MeV;
    traccc::host::kalman_fitting_algorithm host_fitting(fit_cfg, host_mr, copy);

    // Iterate over events
    for (std::size_t i_evt = 0; i_evt < n_events; i_evt++) {

        // Truth Track Candidates
        traccc::event_data evt_data(path, i_evt, host_mr);

        traccc::edm::measurement_collection::host truth_measurements{host_mr};
        traccc::edm::track_container<traccc::default_algebra>::host
            truth_track_candidates{host_mr};
        evt_data.generate_truth_candidates(truth_track_candidates,
                                           truth_measurements, sg, host_mr);
        truth_track_candidates.measurements =
            vecmem::get_data(truth_measurements);

        ASSERT_EQ(truth_track_candidates.tracks.size(), n_truth_tracks);

        // Prepare truth seeds
        traccc::bound_track_parameters_collection_types::host seeds(&host_mr);
        for (unsigned int i_trk = 0; i_trk < n_truth_tracks; i_trk++) {
            seeds.push_back(truth_track_candidates.tracks.at(i_trk).params());
        }
        ASSERT_EQ(seeds.size(), n_truth_tracks);

        // Read measurements
        traccc::edm::measurement_collection::host measurements_per_event{
            host_mr};
        traccc::io::read_measurements(measurements_per_event, i_evt, path);

        // Run finding
        auto track_candidates = host_finding(
            detector, field, vecmem::get_data(measurements_per_event),
            vecmem::get_data(seeds));

        ASSERT_EQ(track_candidates.tracks.size(), n_truth_tracks);

        for (unsigned int i_trk = 0; i_trk < n_truth_tracks; i_trk++) {

            consistency_tests(track_candidates.tracks.at(i_trk),
                              track_candidates.states);

            ndf_tests(track_candidates.tracks.at(i_trk),
                      track_candidates.states, measurements_per_event);
        }

        // Run fitting
        auto track_states = host_fitting(
            detector, field,
            traccc::edm::track_container<traccc::default_algebra>::const_data(
                track_candidates));
        const std::size_t n_fitted_tracks =
            count_successfully_fitted_tracks(track_states.tracks);

        ASSERT_EQ(track_states.tracks.size(), n_truth_tracks);
        ASSERT_EQ(track_states.tracks.size(), n_fitted_tracks);

        for (unsigned int i_trk = 0; i_trk < n_truth_tracks; i_trk++) {

            consistency_tests(track_states.tracks.at(i_trk),
                              track_states.states);

            ndf_tests(track_states.tracks.at(i_trk), track_states.states,
                      measurements_per_event);

            fit_performance_writer.write(
                track_states.tracks.at(i_trk), track_states.states,
                measurements_per_event, detector.as<detector_traits>(),
                evt_data);
        }
    }

    fit_performance_writer.finalize();

    /********************
     * Pull value test
     ********************/

    static const std::vector<std::string> pull_names{
        "pull_d0", "pull_z0", "pull_phi", "pull_theta", "pull_qop"};
    pull_value_tests(fit_writer_cfg.file_path, pull_names);

    /********************
     * P-value test
     ********************/

    p_value_tests(fit_writer_cfg.file_path);

    /********************
     * Success rate test
     ********************/

    float success_rate = static_cast<float>(n_success) /
                         static_cast<float>(n_truth_tracks * n_events);

    ASSERT_FLOAT_EQ(success_rate, 1.00f);
}

INSTANTIATE_TEST_SUITE_P(
    CkfSparseTrackTelescopeValidation0, CkfSparseTrackTelescopeTests,
    ::testing::Values(std::make_tuple(
        "telescope_single_tracks", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 400.f, 400.f},
        std::array<scalar, 2u>{1.f, 1.f}, std::array<scalar, 2u>{0.f, 0.f},
        std::array<scalar, 2u>{0.f, 0.f}, traccc::muon<scalar>(), 1, 5000,
        false, 20.f, 9u, 20.f, vector3{0, 0, 2 * traccc::unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CkfSparseTrackTelescopeValidation1, CkfSparseTrackTelescopeTests,
    ::testing::Values(std::make_tuple(
        "telescope_double_tracks", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 400.f, 400.f},
        std::array<scalar, 2u>{1.f, 1.f}, std::array<scalar, 2u>{0.f, 0.f},
        std::array<scalar, 2u>{0.f, 0.f}, traccc::muon<scalar>(), 2, 2500,
        false, 20.f, 9u, 20.f, vector3{0, 0, 2 * traccc::unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CkfSparseTrackTelescopeValidation2, CkfSparseTrackTelescopeTests,
    ::testing::Values(std::make_tuple(
        "telescope_quadra_tracks", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 400.f, 400.f},
        std::array<scalar, 2u>{1.f, 1.f}, std::array<scalar, 2u>{0.f, 0.f},
        std::array<scalar, 2u>{0.f, 0.f}, traccc::muon<scalar>(), 4, 1250,
        false, 20.f, 9u, 20.f, vector3{0, 0, 2 * traccc::unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CkfSparseTrackTelescopeValidation3, CkfSparseTrackTelescopeTests,
    ::testing::Values(std::make_tuple(
        "telescope_decade_tracks", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 400.f, 400.f},
        std::array<scalar, 2u>{1.f, 1.f}, std::array<scalar, 2u>{0.f, 0.f},
        std::array<scalar, 2u>{0.f, 0.f}, traccc::muon<scalar>(), 10, 500,
        false, 20.f, 9u, 20.f, vector3{0, 0, 2 * traccc::unit<scalar>::T})));

INSTANTIATE_TEST_SUITE_P(
    CkfSparseTrackTelescopeValidation4, CkfSparseTrackTelescopeTests,
    ::testing::Values(std::make_tuple(
        "telescope_decade_tracks_random_charge",
        std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 400.f, 400.f},
        std::array<scalar, 2u>{1.f, 1.f}, std::array<scalar, 2u>{0.f, 0.f},
        std::array<scalar, 2u>{0.f, 0.f}, traccc::muon<scalar>(), 10, 500, true,
        20.f, 9u, 20.f, vector3{0, 0, 2 * traccc::unit<scalar>::T})));
