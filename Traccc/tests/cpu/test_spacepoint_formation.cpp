/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "tests/test_detectors.hpp"
#include "traccc/definitions/common.hpp"
#include "traccc/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace traccc;

TEST(spacepoint_formation, cpu) {

    // Memory resource used by the EDM.
    vecmem::host_memory_resource host_mr;

    // Use rectangle surfaces
    detray::mask<detray::rectangle2D, traccc::default_algebra> rectangle{
        0u, 10000.f * traccc::unit<scalar>::mm,
        10000.f * traccc::unit<scalar>::mm};

    // Plane alignment direction (aligned to x-axis)
    detray::detail::ray<traccc::default_algebra> traj{
        {0, 0, 0}, 0, {1, 0, 0}, -1};
    // Position of planes (in mm unit)
    std::vector<scalar> plane_positions = {20.f,  40.f,  60.f,  80.f, 100.f,
                                           120.f, 140.f, 160.f, 180.f};

    detray::tel_det_config tel_cfg{rectangle};
    tel_cfg.positions(plane_positions);
    tel_cfg.pilot_track(traj);

    // Create telescope geometry
    auto [det, name_map] = build_telescope_detector(host_mr, tel_cfg);

    auto surfaces = det.surfaces();

    traccc::host_detector host_det;
    host_det.set<traccc::telescope_detector>(std::move(det));

    // Surface lookup

    // Prepare measurement collection
    edm::measurement_collection::host measurements{host_mr};

    // Add a measurement at the first plane
    measurements.push_back({{7.f, 2.f},
                            {0.f, 0.f},
                            2,
                            0.f,
                            0.f,
                            0u,
                            surfaces[0].identifier(),
                            {1u, 1u},
                            0u});

    // Add a measurement at the last plane
    measurements.push_back({{10.f, 15.f},
                            {0.f, 0.f},
                            2u,
                            0.f,
                            0.f,
                            0u,
                            surfaces[8u].identifier(),
                            {1u, 1u},
                            1u});

    // Run spacepoint formation
    host::silicon_pixel_spacepoint_formation_algorithm sp_formation(host_mr);
    auto spacepoints = sp_formation(host_det, vecmem::get_data(measurements));

    // Check the results
    EXPECT_EQ(spacepoints.size(), 2u);
    EXPECT_FLOAT_EQ(static_cast<float>(spacepoints[0].x()), 20.f);
    EXPECT_FLOAT_EQ(static_cast<float>(spacepoints[0].y()), 7.f);
    EXPECT_FLOAT_EQ(static_cast<float>(spacepoints[0].z()), 2.f);
    EXPECT_FLOAT_EQ(static_cast<float>(spacepoints[1].x()), 180.f);
    EXPECT_FLOAT_EQ(static_cast<float>(spacepoints[1].y()), 10.f);
    EXPECT_FLOAT_EQ(static_cast<float>(spacepoints[1].z()), 15.f);
}
