/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "tests/test_detectors.hpp"
#include "traccc/definitions/common.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/sycl/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/sycl/shared_memory_resource.hpp>
#include <vecmem/utils/sycl/copy.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace traccc;

TEST(SYCLSpacepointFormation, sycl) {

    // Creating SYCL queue object
    vecmem::sycl::queue_wrapper vecmem_queue;
    traccc::sycl::queue_wrapper traccc_queue(vecmem_queue.queue());
    std::cout << "Running on device: " << vecmem_queue.device_name() << "\n";

    // Memory resource used by the EDM.
    vecmem::sycl::shared_memory_resource shared_mr{vecmem_queue};
    traccc::memory_resource mr{shared_mr};

    // Copy object
    vecmem::sycl::copy copy{vecmem_queue};
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
    auto [det, name_map] = build_telescope_detector(shared_mr, tel_cfg);

    // Surface lookup
    auto surfaces = det.surfaces();

    traccc::host_detector polymorphic_detector;
    polymorphic_detector.set<traccc::telescope_detector>(std::move(det));

    const traccc::detector_buffer detector_buffer =
        traccc::buffer_from_host_detector(polymorphic_detector, shared_mr,
                                          copy);

    // Prepare measurement collection
    edm::measurement_collection::host measurements{shared_mr};

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
    traccc::sycl::silicon_pixel_spacepoint_formation_algorithm sp_formation(
        mr, copy, traccc_queue);
    auto spacepoints_buffer =
        sp_formation(detector_buffer, vecmem::get_data(measurements));

    edm::spacepoint_collection::device spacepoints(spacepoints_buffer);

    // Check the results
    EXPECT_EQ(copy.get_size(spacepoints_buffer), 2u);
    std::set<point3> test;
    test.insert(spacepoints[0].global());
    test.insert(spacepoints[1].global());

    std::set<point3> ref;
    ref.insert({180.f, 10.f, 15.f});
    ref.insert({20.f, 7.f, 2.f});

    EXPECT_EQ(test, ref);
}
