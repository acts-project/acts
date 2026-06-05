/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "tests/test_detectors.hpp"
#include "traccc/bfield/construct_const_bfield.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/io/csv/make_hit_reader.hpp"
#include "traccc/io/csv/make_measurement_hit_id_reader.hpp"
#include "traccc/io/csv/make_measurement_reader.hpp"
#include "traccc/io/csv/make_particle_reader.hpp"
#include "traccc/simulation/event_generators.hpp"
#include "traccc/simulation/simulator.hpp"

// Detray include(s).
#include <detray/geometry/mask.hpp>
#include <detray/geometry/shapes/line.hpp>
#include <detray/geometry/shapes/rectangle2D.hpp>
#include <detray/geometry/tracking_surface.hpp>
#include <detray/test/utils/statistics.hpp>

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <filesystem>

using namespace traccc;

constexpr scalar tol{1e-7f};

TEST(traccc_simulation, simulation) {

    using line_t = detray::mask<detray::line<false>, traccc::default_algebra>;
    using rectangle_t =
        detray::mask<detray::rectangle2D, traccc::default_algebra>;

    traccc::bound_track_parameters<traccc::default_algebra> bound_params{};
    bound_params.set_bound_local({1.f, 2.f});

    measurement_smearer<traccc::default_algebra> smearer(0.f, 0.f);

    traccc::io::csv::measurement iomeas1;
    smearer.template operator()<line_t>({-3.f, 2.f}, bound_params, iomeas1);
    ASSERT_NEAR(iomeas1.local0, 0.f, tol);
    ASSERT_NEAR(iomeas1.local1, 0.f, tol);

    traccc::io::csv::measurement iomeas2;
    smearer.template operator()<line_t>({2.f, -5.f}, bound_params, iomeas2);
    ASSERT_NEAR(iomeas2.local0, 3.f, tol);
    ASSERT_NEAR(iomeas2.local1, 0.f, tol);

    traccc::io::csv::measurement iomeas3;
    smearer.template operator()<rectangle_t>({2.f, -5.f}, bound_params,
                                             iomeas3);
    ASSERT_NEAR(iomeas3.local0, 3.f, tol);
    ASSERT_NEAR(iomeas3.local1, -3.f, tol);
}

GTEST_TEST(traccc_simulation, toy_detector_simulation) {

    // Create geometry
    vecmem::host_memory_resource host_mr;

    // Create B field
    using b_field_t = covfie::field<traccc::const_bfield_backend_t<scalar>>;
    const vector3 B{0.f, 0.f, 2.f * traccc::unit<scalar>::T};
    b_field_t field =
        traccc::construct_const_bfield(B)
            .as_field<traccc::const_bfield_backend_t<traccc::scalar>>();

    // Create geometry
    detray::toy_det_config<scalar> toy_cfg{};
    const auto [detector, names] =
        detray::build_toy_detector<traccc::default_algebra>(host_mr, toy_cfg);

    using geo_cxt_t = typename decltype(detector)::geometry_context;
    const geo_cxt_t ctx{};

    // Create track generator
    using uniform_gen_t =
        detray::detail::random_numbers<scalar,
                                       std::uniform_real_distribution<scalar>>;
    using generator_type =
        detray::random_track_generator<traccc::free_track_parameters<>,
                                       uniform_gen_t>;
    generator_type::configuration gen_cfg{};
    constexpr unsigned int n_tracks{2500u};
    const vector3 ori{0.f, 0.f, 0.f};
    gen_cfg.n_tracks(n_tracks);
    gen_cfg.origin(ori);
    // @TODO The simulator sometimes gets stuck for lower momentum
    gen_cfg.p_tot(5.f * traccc::unit<scalar>::GeV);
    generator_type generator(gen_cfg);

    // Create smearer
    measurement_smearer<traccc::default_algebra> smearer(
        67.f * traccc::unit<scalar>::um, 170.f * traccc::unit<scalar>::um);

    std::size_t n_events{10u};

    using detector_type = decltype(detector);
    using writer_type =
        smearing_writer<measurement_smearer<traccc::default_algebra>>;

    typename writer_type::config writer_cfg{smearer};

    auto sim = simulator<detector_type, b_field_t, generator_type, writer_type>(
        traccc::muon<scalar>(), n_events, detector, field, std::move(generator),
        std::move(writer_cfg));

    // Lift step size constraints
    sim.get_config().propagation.stepping.step_constraint =
        std::numeric_limits<float>::max();
    sim.get_config().propagation.navigation.search_window = {3u, 3u};

    // Do the simulation
    sim.run();

    for (std::size_t i_event = 0u; i_event < n_events; i_event++) {

        std::vector<traccc::io::csv::particle> particles;
        auto particle_reader = traccc::io::csv::make_particle_reader(
            traccc::io::get_event_filename(i_event, "-particles_initial.csv"));
        traccc::io::csv::particle io_particle;
        while (particle_reader.read(io_particle)) {
            particles.push_back(io_particle);
        }

        std::vector<traccc::io::csv::hit> hits;
        auto hit_reader = traccc::io::csv::make_hit_reader(
            traccc::io::get_event_filename(i_event, "-hits.csv"));
        traccc::io::csv::hit io_hit;
        while (hit_reader.read(io_hit)) {
            hits.push_back(io_hit);
        }

        std::vector<traccc::io::csv::measurement> measurements;
        auto measurement_reader = traccc::io::csv::make_measurement_reader(
            traccc::io::get_event_filename(i_event, "-measurements.csv"));
        traccc::io::csv::measurement io_measurement;
        while (measurement_reader.read(io_measurement)) {
            measurements.push_back(io_measurement);
        }

        std::vector<traccc::io::csv::measurement_hit_id> meas_hit_ids;
        auto measurement_hit_id_reader =
            traccc::io::csv::make_measurement_hit_id_reader(
                traccc::io::get_event_filename(i_event,
                                               "-measurement-simhit-map.csv"));
        traccc::io::csv::measurement_hit_id io_meas_hit_id;
        while (measurement_hit_id_reader.read(io_meas_hit_id)) {
            meas_hit_ids.push_back(io_meas_hit_id);
        }

        ASSERT_EQ(particles.size(), n_tracks);
        ASSERT_TRUE(not measurements.empty());
        ASSERT_EQ(hits.size(), measurements.size());
        ASSERT_EQ(hits.size(), meas_hit_ids.size());

        // Let's check if measurement smearing works correctly...
        std::vector<scalar> local0_diff;
        std::vector<scalar> local1_diff;

        const std::size_t nhits = hits.size();
        for (std::size_t i = 0u; i < nhits; i++) {
            const point3 pos{hits[i].tx, hits[i].ty, hits[i].tz};
            const vector3 mom{hits[i].tpx, hits[i].tpy, hits[i].tpz};
            const auto truth_local =
                detray::tracking_surface{
                    detector, detray::geometry::identifier(hits[i].geometry_id)}
                    .global_to_local(ctx, pos, vector::normalize(mom));

            local0_diff.push_back(truth_local[0] - measurements[i].local0);
            local1_diff.push_back(truth_local[1] - measurements[i].local1);

            ASSERT_EQ(meas_hit_ids[i].hit_id, i);
            ASSERT_EQ(meas_hit_ids[i].measurement_id, i);
        }

        const auto var0 = detray::statistics::variance(local0_diff);
        const auto var1 = detray::statistics::variance(local1_diff);

        EXPECT_NEAR((std::sqrt(var0) - smearer.stddev[0]) / smearer.stddev[0],
                    0.f, 0.1f);
        EXPECT_NEAR((std::sqrt(var1) - smearer.stddev[1]) / smearer.stddev[1],
                    0.f, 0.1f);
    }
}

// Test parameters: <initial momentum, theta direction, charge>
class TelescopeDetectorSimulation
    : public ::testing::TestWithParam<
          std::tuple<std::string, scalar, scalar, scalar>> {};

TEST_P(TelescopeDetectorSimulation, telescope_detector_simulation) {

    // Create geometry
    vecmem::host_memory_resource host_mr;

    // Build from given module positions
    std::vector<scalar> positions = {0.f,   50.f,  100.f, 150.f, 200.f, 250.f,
                                     300.f, 350.f, 400.f, 450.f, 500.f};

    // A thickness larger than 0.1 cm will flip the track direction of low
    // energy (or non-relativistic) particle due to the large scattering
    const scalar thickness = 0.005f * traccc::unit<scalar>::cm;

    detray::tel_det_config<traccc::default_algebra, detray::rectangle2D>
        tel_cfg{1000.f * traccc::unit<scalar>::mm,
                1000.f * traccc::unit<scalar>::mm};
    tel_cfg.positions(positions).mat_thickness(thickness);

    const auto [detector, names] =
        detray::build_telescope_detector(host_mr, tel_cfg);

    // Directory name
    const std::string directory = std::get<0>(GetParam()) + "/";
    std::filesystem::create_directory(directory);

    // Field
    using b_field_t = covfie::field<traccc::const_bfield_backend_t<scalar>>;
    const vector3 B{0.f, 0.f, 2.f * traccc::unit<scalar>::T};
    b_field_t field =
        traccc::construct_const_bfield(B)
            .as_field<traccc::const_bfield_backend_t<traccc::scalar>>();

    // Momentum
    const scalar mom = std::get<1>(GetParam());

    // Create track generator
    constexpr unsigned int theta_steps{1u};
    constexpr unsigned int phi_steps{1u};
    const vector3 ori{0.f, 0.f, 0.f};
    const scalar theta = std::get<2>(GetParam());

    const scalar charge = std::get<3>(GetParam());

    // Track generator
    using generator_type =
        detray::uniform_track_generator<traccc::free_track_parameters<>>;
    generator_type::configuration gen_cfg{};
    gen_cfg.theta_steps(theta_steps);
    gen_cfg.phi_steps(phi_steps);
    gen_cfg.origin(ori);
    gen_cfg.theta_range(theta, theta);
    gen_cfg.p_tot(mom);
    gen_cfg.charge(charge);
    generator_type generator(gen_cfg);

    // Create smearer
    measurement_smearer<traccc::default_algebra> smearer(
        50.f * traccc::unit<scalar>::um, 50.f * traccc::unit<scalar>::um);

    std::size_t n_events{1000u};

    using detector_type = decltype(detector);
    using generator_type = decltype(generator);
    using writer_type =
        smearing_writer<measurement_smearer<traccc::default_algebra>>;

    typename writer_type::config writer_cfg{smearer};

    auto sim = simulator<detector_type, b_field_t, generator_type, writer_type>(
        traccc::muon<scalar>(), n_events, detector, field, std::move(generator),
        std::move(writer_cfg), directory);

    // Lift step size constraints
    sim.get_config().propagation.stepping.step_constraint =
        std::numeric_limits<float>::max();

    // Run simulation
    sim.get_config().propagation.navigation.intersection.overstep_tolerance =
        -100.f * unit<float>::um;
    sim.get_config().propagation.navigation.intersection.max_mask_tolerance =
        1.f * unit<float>::mm;
    sim.run();

    for (std::size_t i_event{0u}; i_event < n_events; i_event++) {

        std::vector<traccc::io::csv::measurement> measurements;
        auto measurement_reader = traccc::io::csv::make_measurement_reader(
            directory +
            traccc::io::get_event_filename(i_event, "-measurements.csv"));
        traccc::io::csv::measurement io_measurement;
        while (measurement_reader.read(io_measurement)) {
            measurements.push_back(io_measurement);
        }

        // Make sure that number of measurements is equal to the number of
        // physical planes
        ASSERT_EQ(measurements.size(), positions.size());
    }
}

INSTANTIATE_TEST_SUITE_P(
    Simulation, TelescopeDetectorSimulation,
    ::testing::Values(
        std::make_tuple("0", 0.1f * traccc::unit<scalar>::GeV, 0.01f, -1.f),
        std::make_tuple("1", 1.f * traccc::unit<scalar>::GeV, 0.01f, -1.f),
        std::make_tuple("2", 10.f * traccc::unit<scalar>::GeV, 0.01f, -1.f),
        std::make_tuple("3", 100.f * traccc::unit<scalar>::GeV, 0.01f, -1.f),
        std::make_tuple("4", 0.1f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 12.f, 1.f),
        std::make_tuple("5", 1.f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 12.f, 1.f),
        std::make_tuple("6", 10.f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 12.f, 1.f),
        std::make_tuple("7", 100.f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 12.f, 1.f),
        std::make_tuple("8", 0.1f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 8.f, -1.f),
        std::make_tuple("9", 1.f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 8.f, -1.f),
        std::make_tuple("10", 10.f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 8.f, -1.f),
        std::make_tuple("11", 100.f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 8.f, -1.f),
        std::make_tuple("12", 0.1f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 6.f, 1.f),
        std::make_tuple("13", 1.f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 6.f, 1.f),
        std::make_tuple("14", 10.f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 6.f, 1.f),
        std::make_tuple("15", 100.f * traccc::unit<scalar>::GeV,
                        traccc::constant<scalar>::pi / 6.f, 1.f)));
