/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/bfield/construct_const_bfield.hpp"
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
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

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <filesystem>
#include <string>

using namespace traccc;
// This defines the local frame test suite
TEST_P(CpuCkfCombinatoricsTelescopeTests, Run) {

    // Get the parameters
    const std::string name = std::get<0>(GetParam());
    const std::array<scalar, 3u> origin = std::get<1>(GetParam());
    const std::array<scalar, 3u> origin_stddev = std::get<2>(GetParam());
    const std::array<scalar, 2u> mom_range = std::get<3>(GetParam());
    const std::array<scalar, 2u> eta_range = std::get<4>(GetParam());
    const std::array<scalar, 2u> theta_range = eta_to_theta_range(eta_range);
    const std::array<scalar, 2u> phi_range = std::get<5>(GetParam());
    const unsigned int n_truth_tracks = std::get<7>(GetParam());
    const unsigned int n_events = std::get<8>(GetParam());
    const bool random_charge = std::get<9>(GetParam());

    /*****************************
     * Build a telescope geometry
     *****************************/

    // Memory resources used by the application.
    vecmem::host_memory_resource host_mr;

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

    const auto field = traccc::construct_const_bfield(std::get<13>(GetParam()));

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
        std::get<6>(GetParam()), n_events, detector.as<detector_traits>(),
        field.as_field<const_bfield_backend_t<traccc::scalar>>(),
        std::move(generator), std::move(smearer_writer_cfg), full_path);
    sim.run();

    /*****************************
     * Do the reconstruction
     *****************************/

    // Seed generator
    seed_generator<host_detector_type> sg(detector.as<detector_traits>(),
                                          stddevs);

    // Finding algorithm configuration
    traccc::finding_config cfg_no_limit;
    cfg_no_limit.max_num_branches_per_seed =
        std::numeric_limits<unsigned int>::max();
    cfg_no_limit.max_num_branches_per_surface = 10;
    cfg_no_limit.chi2_max = 30.f;
    cfg_no_limit.duplicate_removal_minimum_length = 100u;

    traccc::finding_config cfg_limit;
    cfg_limit.max_num_branches_per_seed = 500;
    cfg_limit.max_num_branches_per_surface = 10;
    cfg_limit.chi2_max = 30.f;
    cfg_limit.duplicate_removal_minimum_length = 100u;

    // Finding algorithm object
    traccc::host::combinatorial_kalman_filter_algorithm host_finding(
        cfg_no_limit, host_mr);
    traccc::host::combinatorial_kalman_filter_algorithm host_finding_limit(
        cfg_limit, host_mr);

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
        const traccc::bound_track_parameters_collection_types::const_view
            seeds_view = vecmem::get_data(seeds);

        // Read measurements
        traccc::edm::measurement_collection::host measurements_per_event{
            host_mr};
        traccc::io::read_measurements(measurements_per_event, i_evt, path);
        const traccc::edm::measurement_collection::const_view
            measurements_view = vecmem::get_data(measurements_per_event);

        // Run finding
        auto track_candidates =
            host_finding(detector, field, measurements_view, seeds_view);

        auto track_candidates_limit =
            host_finding_limit(detector, field, measurements_view, seeds_view);

        // Make sure that the number of found tracks = n_track ^ (n_planes + 1)
        ASSERT_TRUE(track_candidates.tracks.size() >
                    track_candidates_limit.tracks.size());
        ASSERT_EQ(track_candidates.tracks.size(),
                  std::pow(n_truth_tracks, std::get<11>(GetParam()) + 1));
        ASSERT_EQ(track_candidates_limit.tracks.size(),
                  n_truth_tracks * cfg_limit.max_num_branches_per_seed);
    }
}

// Testing two identical tracks
INSTANTIATE_TEST_SUITE_P(
    CpuCkfCombinatoricsTelescopeValidation0, CpuCkfCombinatoricsTelescopeTests,
    ::testing::Values(std::make_tuple(
        "telescope_combinatorics_twin", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 2u>{100.f, 100.f}, std::array<scalar, 2u>{0.f, 0.f},
        std::array<scalar, 2u>{0.f, 0.f}, traccc::muon<scalar>(), 2, 1, false,
        20.f, 9u, 20.f, vector3{2 * traccc::unit<scalar>::T, 0, 0})));

// Testing three identical tracks
INSTANTIATE_TEST_SUITE_P(
    CpuCkfCombinatoricsTelescopeValidation1, CpuCkfCombinatoricsTelescopeTests,
    ::testing::Values(std::make_tuple(
        "telescope_combinatorics_trio", std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 3u>{0.f, 0.f, 0.f},
        std::array<scalar, 2u>{100.f, 100.f}, std::array<scalar, 2u>{0.f, 0.f},
        std::array<scalar, 2u>{0.f, 0.f}, traccc::muon<scalar>(), 3, 1, false,
        20.f, 9u, 20.f, vector3{2 * traccc::unit<scalar>::T, 0, 0})));
