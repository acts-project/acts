/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/bfield/construct_const_bfield.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/resolution/fitting_performance_writer.hpp"
#include "traccc/simulation/event_generators.hpp"
#include "traccc/simulation/measurement_smearer.hpp"
#include "traccc/simulation/simulator.hpp"
#include "traccc/simulation/smearing_writer.hpp"
#include "traccc/utils/ranges.hpp"
#include "traccc/utils/seed_generator.hpp"

// Test include(s).
#include "tests/kalman_fitting_wire_chamber_test.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <filesystem>
#include <string>

using namespace traccc;

TEST_P(KalmanFittingWireChamberTests, Run) {

    // Get the parameters
    const std::string name = std::get<0>(GetParam());
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
     * Build a drift chamber
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
            std::filesystem::path(path + "wire_chamber_geometry.json"))
            .native(),
        std::filesystem::absolute(
            std::filesystem::path(path +
                                  "wire_chamber_homogeneous_material.json"))
            .native(),
        std::filesystem::absolute(
            std::filesystem::path(path + "wire_chamber_surface_grids.json"))
            .native());
    const auto field = traccc::construct_const_bfield(B);

    /***************************
     * Generate simulation data
     ***************************/

    // Track generator
    using generator_type =
        detray::random_track_generator<traccc::free_track_parameters<>,
                                       uniform_gen_t>;
    generator_type::configuration gen_cfg{};
    gen_cfg.n_tracks(n_truth_tracks);
    gen_cfg.origin(std::get<1>(GetParam()));
    gen_cfg.origin_stddev(std::get<2>(GetParam()));
    gen_cfg.phi_range(std::get<5>(GetParam()));
    gen_cfg.eta_range(std::get<4>(GetParam()));
    gen_cfg.mom_range(std::get<3>(GetParam()));
    gen_cfg.randomize_charge(random_charge);
    gen_cfg.seed(42);
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

    sim.get_config().propagation.navigation.search_window = search_window;

    sim.run();

    /***************
     * Run fitting
     ***************/

    // Seed generator
    seed_generator<host_detector_type> sg(detector.as<detector_traits>(),
                                          stddevs);

    // Fitting algorithm object
    traccc::fitting_config fit_cfg;
    fit_cfg.propagation.navigation.intersection.min_mask_tolerance =
        static_cast<float>(mask_tolerance);
    fit_cfg.propagation.navigation.search_window = search_window;
    // TODO: Disable until overlaps are handled correctly
    fit_cfg.propagation.navigation.estimate_scattering_noise = false;
    fit_cfg.ptc_hypothesis = ptc;
    fit_cfg.min_pT = 100.f * traccc::unit<float>::MeV;
    traccc::host::kalman_fitting_algorithm fitting(fit_cfg, host_mr, copy);

    // Iterate over events
    for (std::size_t i_evt = 0; i_evt < n_events; i_evt++) {

        // Event map
        traccc::event_data evt_data(path, i_evt, host_mr);
        // Truth Track Candidates
        traccc::edm::measurement_collection::host measurements(host_mr);
        traccc::edm::track_container<traccc::default_algebra>::host
            track_candidates{host_mr};
        evt_data.generate_truth_candidates(track_candidates, measurements, sg,
                                           host_mr);
        track_candidates.measurements = vecmem::get_data(measurements);

        // n_trakcs = 100
        ASSERT_EQ(track_candidates.tracks.size(), n_truth_tracks);

        // Run fitting
        auto track_states = fitting(
            detector, field,
            traccc::edm::track_container<traccc::default_algebra>::const_data(
                track_candidates));

        // Iterator over tracks
        const std::size_t n_tracks = track_states.tracks.size();

        ASSERT_GE(static_cast<float>(n_tracks),
                  0.98 * static_cast<float>(n_truth_tracks));

        const std::size_t n_fitted_tracks =
            count_successfully_fitted_tracks(track_states.tracks);
        ASSERT_GE(static_cast<float>(n_fitted_tracks),
                  0.93f * static_cast<float>(n_truth_tracks));

        for (std::size_t i_trk = 0; i_trk < n_tracks; i_trk++) {

            // Some fits fail. The results of those cannot be reasonably tested.
            if (track_states.tracks.at(i_trk).fit_outcome() !=
                traccc::track_fit_outcome::SUCCESS) {
                continue;
            }

            consistency_tests(track_states.tracks.at(i_trk),
                              track_states.states);

            ndf_tests(track_states.tracks.at(i_trk), track_states.states,
                      measurements);

            fit_performance_writer.write(
                track_states.tracks.at(i_trk), track_states.states,
                measurements, detector.as<detector_traits>(), evt_data);
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

    //@TODO: Develop an extension of KF-based fitter (e.g. Deterministic
    // Annealing Filter) to resolve left-right ambiguity and pass the p-value
    // test
    // p_value_tests(fit_writer_cfg.file_path);

    /********************
     * Success rate test
     ********************/

    scalar success_rate = static_cast<scalar>(n_success) /
                          static_cast<scalar>(n_truth_tracks * n_events);

    // TODO: Raise back to 95%
    ASSERT_GE(success_rate, 0.93f);
    ASSERT_LE(success_rate, 1.00f);
}

INSTANTIATE_TEST_SUITE_P(
    KalmanFitWireChamberValidation0, KalmanFittingWireChamberTests,
    ::testing::Values(std::make_tuple(
        "wire_2_GeV_muon", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f}, std::array<scalar, 2u>{2.f, 2.f},
        std::array<scalar, 2u>{-1.f, 1.f},
        std::array<scalar, 2u>{-traccc::constant<scalar>::pi,
                               traccc::constant<scalar>::pi},
        traccc::muon<scalar>(), 100, 100, false)));

// @TODO: Make full eta range work
INSTANTIATE_TEST_SUITE_P(
    KalmanFitWireChamberValidation1, KalmanFittingWireChamberTests,
    ::testing::Values(std::make_tuple(
        "wire_10_GeV_muon", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 2u>{10.f, 10.f}, std::array<scalar, 2u>{-0.3f, 0.3f},
        std::array<scalar, 2u>{-traccc::constant<scalar>::pi,
                               traccc::constant<scalar>::pi},
        traccc::muon<scalar>(), 100, 100, false)));

// @TODO: Make full eta range work
INSTANTIATE_TEST_SUITE_P(
    KalmanFitWireChamberValidation2, KalmanFittingWireChamberTests,
    ::testing::Values(std::make_tuple(
        "wire_100_GeV_muon", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 2u>{100.f, 100.f},
        std::array<scalar, 2u>{-0.4f, 0.4f},
        std::array<scalar, 2u>{-traccc::constant<scalar>::pi,
                               traccc::constant<scalar>::pi},
        traccc::muon<scalar>(), 100, 100, false)));

INSTANTIATE_TEST_SUITE_P(
    KalmanFitWireChamberValidation3, KalmanFittingWireChamberTests,
    ::testing::Values(std::make_tuple(
        "wire_2_GeV_anti_muon", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f}, std::array<scalar, 2u>{2.f, 2.f},
        std::array<scalar, 2u>{-1.f, 1.f},
        std::array<scalar, 2u>{-traccc::constant<scalar>::pi,
                               traccc::constant<scalar>::pi},
        traccc::antimuon<scalar>(), 100, 100, false)));

INSTANTIATE_TEST_SUITE_P(
    KalmanFitWireChamberValidation4, KalmanFittingWireChamberTests,
    ::testing::Values(std::make_tuple(
        "wire_2_GeV_random_charge", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f}, std::array<scalar, 2u>{2.f, 2.f},
        std::array<scalar, 2u>{-1.f, 1.f},
        std::array<scalar, 2u>{-traccc::constant<scalar>::pi,
                               traccc::constant<scalar>::pi},
        traccc::antimuon<scalar>(), 100, 100, true)));
