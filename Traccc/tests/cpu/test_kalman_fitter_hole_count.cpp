/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
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
#include "traccc/simulation/simulator.hpp"
#include "traccc/utils/ranges.hpp"
#include "traccc/utils/seed_generator.hpp"

// Test include(s).
#include "tests/kalman_fitting_telescope_test.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <filesystem>
#include <string>

using namespace traccc;

class KalmanFittingHoleCountTests : public KalmanFittingTelescopeTests {};

TEST_P(KalmanFittingHoleCountTests, Run) {

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

    // We only test one track of one event
    ASSERT_EQ(n_truth_tracks, 1u);
    ASSERT_EQ(n_events, 1u);

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

    /***************
     * Run fitting
     ***************/

    // Seed generator
    seed_generator<host_detector_type> sg(detector.as<detector_traits>(),
                                          stddevs);

    // Fitting algorithm object
    traccc::fitting_config fit_cfg;
    fit_cfg.ptc_hypothesis = ptc;
    fit_cfg.min_p = 10.f * traccc::unit<float>::MeV;
    fit_cfg.min_pT = 60.f * traccc::unit<float>::MeV;
    traccc::host::kalman_fitting_algorithm fitting(fit_cfg, host_mr, copy);

    // Event map
    traccc::event_data evt_data(path, 0u, host_mr);

    // Truth Track Candidates
    traccc::edm::measurement_collection::host measurements(host_mr);
    traccc::edm::track_container<traccc::default_algebra>::host
        track_candidates{host_mr};
    evt_data.generate_truth_candidates(track_candidates, measurements, sg,
                                       host_mr);
    track_candidates.measurements = vecmem::get_data(measurements);
    // Measurement index vector
    auto& cands = track_candidates.tracks.at(0u).constituent_links();

    // Some sanity checks
    ASSERT_EQ(track_candidates.tracks.size(), n_truth_tracks);
    const auto n_planes = std::get<11>(GetParam());
    ASSERT_EQ(cands.size(), n_planes);

    // Pop some track candidates to create holes
    // => The number of holes = 8
    ASSERT_TRUE(cands.size() > 8u);
    cands.erase(cands.begin());
    cands.erase(cands.begin());
    cands.erase(cands.begin() + 2);
    cands.erase(cands.begin() + 2);
    cands.erase(cands.begin() + 7);
    cands.pop_back();
    cands.pop_back();
    cands.pop_back();

    // A sanity check on the number of candidiates
    ASSERT_EQ(cands.size(), n_planes - 8u);

    // Run fitting
    auto track_states = fitting(
        detector, field,
        traccc::edm::track_container<traccc::default_algebra>::const_data(
            track_candidates));

    // A sanity check
    const std::size_t n_tracks = track_states.tracks.size();
    ASSERT_EQ(n_tracks, n_truth_tracks);

    // Check the number of holes
    // The three holes at the end are not counted as KF aborts once it goes
    // through all track candidates
    const auto track = track_states.tracks.at(0u);
    ASSERT_EQ(track.nholes(), 5u);

    // Some sanity checks
    ASSERT_FLOAT_EQ(
        static_cast<float>(track.ndf()),
        static_cast<float>(track.constituent_links().size()) * 2.f - 5.f);
}

INSTANTIATE_TEST_SUITE_P(
    KalmanFittingHoleCount, KalmanFittingHoleCountTests,
    ::testing::Values(std::make_tuple(
        "telescope_1_GeV_0_phi_muon", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f}, std::array<scalar, 2u>{1.f, 1.f},
        std::array<scalar, 2u>{0.f, 0.f}, std::array<scalar, 2u>{0.f, 0.f},
        traccc::muon<scalar>(), 1, 1, false, 20.f, 20u, 20.f,
        vector3{2 * traccc::unit<scalar>::T, 0, 0})));
