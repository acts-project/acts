/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace traccc;

namespace {

// Memory resource used by the EDM.
vecmem::host_memory_resource host_mr;

// Set B field
static constexpr vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                           2.f * unit<scalar>::T};

}  // namespace

// Seeding with two muons
TEST(seeding, case1) {

    // Config objects
    traccc::seedfinder_config finder_config;
    traccc::spacepoint_grid_config grid_config(finder_config);
    traccc::seedfilter_config filter_config;

    // Adjust parameters
    finder_config.deltaRMax = 100.f * unit<float>::mm;
    finder_config.maxPtScattering = 0.5f * unit<float>::GeV;
    traccc::host::seeding_algorithm sa(finder_config, grid_config,
                                       filter_config, host_mr);

    edm::measurement_collection::host measurements(host_mr);
    edm::spacepoint_collection::host spacepoints{host_mr};

    // Spacepoints from 16.62 GeV muon
    measurements.resize(5);
    spacepoints.reserve(5);
    spacepoints.push_back(
        {0,
         traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         {36.6706f, 10.6472f, 104.131f},
         0.f,
         0.f});
    spacepoints.push_back(
        {1,
         traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         {94.2191f, 29.6699f, 113.628f},
         0.f,
         0.f});
    spacepoints.push_back(
        {2,
         traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         {149.805f, 47.9518f, 122.979f},
         0.f,
         0.f});
    spacepoints.push_back(
        {3,
         traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         {218.514f, 70.3049f, 134.029f},
         0.f,
         0.f});
    spacepoints.push_back(
        {4,
         traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         {275.359f, 88.668f, 143.378f},
         0.f,
         0.f});

    // Run seeding
    auto seeds = sa(vecmem::get_data(spacepoints));

    // The number of seeds should be eqaul to one
    ASSERT_EQ(seeds.size(), 1u);

    traccc::track_params_estimation_config track_params_estimation_config;
    traccc::host::track_params_estimation tp(track_params_estimation_config,
                                             host_mr);

    auto bound_params =
        tp(vecmem::get_data(measurements), vecmem::get_data(spacepoints),
           vecmem::get_data(seeds), B);

    // The number of bound track parameters should be eqaul to one
    ASSERT_EQ(bound_params.size(), 1u);

    // Make sure that we have reasonable estimation on momentum
    /* Currently disabled
    EXPECT_NEAR(bound_params[0].p(), 16.62 * unit<scalar>::GeV,
                0.1 * unit<scalar>::GeV);
    */
}

TEST(seeding, case2) {

    // Config objects
    traccc::seedfinder_config finder_config;
    traccc::spacepoint_grid_config grid_config(finder_config);
    traccc::seedfilter_config filter_config;

    // Adjust parameters
    finder_config.deltaRMax = 100.f * unit<float>::mm;
    finder_config.maxPtScattering = 0.5f * unit<float>::GeV;
    traccc::host::seeding_algorithm sa(finder_config, grid_config,
                                       filter_config, host_mr);

    edm::measurement_collection::host measurements(host_mr);
    edm::spacepoint_collection::host spacepoints{host_mr};

    // Spacepoints from 1.85 GeV muon
    measurements.resize(5);
    spacepoints.reserve(5);
    spacepoints.push_back(
        {0,
         traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         {36.301f, 13.1197f, 106.83f},
         0.f,
         0.f});
    spacepoints.push_back(
        {1,
         traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         {93.9366f, 33.7101f, 120.978f},
         0.f,
         0.f});
    spacepoints.push_back(
        {2,
         traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         {149.192f, 52.0562f, 134.678f},
         0.f,
         0.f});
    spacepoints.push_back(
        {3,
         traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         {218.398f, 73.1025f, 151.979f},
         0.f,
         0.f});
    spacepoints.push_back(
        {4,
         traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         {275.322f, 89.0663f, 166.229f},
         0.f,
         0.f});

    // Run seeding
    auto seeds = sa(vecmem::get_data(spacepoints));

    // The number of seeds should be eqaul to one
    ASSERT_EQ(seeds.size(), 1u);

    traccc::track_params_estimation_config track_params_estimation_config;
    traccc::host::track_params_estimation tp(track_params_estimation_config,
                                             host_mr);

    auto bound_params =
        tp(vecmem::get_data(measurements), vecmem::get_data(spacepoints),
           vecmem::get_data(seeds), B);

    // The number of bound track parameters should be eqaul to one
    ASSERT_EQ(bound_params.size(), 1u);

    // Make sure that we have reasonable estimation on momentum
    /* Currently disabled
    EXPECT_NEAR(bound_params[0].p(), 1.85 * unit<scalar>::GeV,
                0.1 * unit<scalar>::GeV);
    */
}
