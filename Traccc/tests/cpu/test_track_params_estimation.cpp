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
#include "traccc/utils/detray_conversion.hpp"

// Detray include(s).
#include <detray/tracks/helix.hpp>

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

using namespace traccc;

namespace {

// Memory resource used by the EDM.
vecmem::host_memory_resource host_mr;

}  // namespace

TEST(track_params_estimation, helix_negative_charge) {

    // Set B field
    const vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                    2.f * unit<scalar>::T};

    // Track property
    const scalar q{-1.f * unit<scalar>::e};
    const point3 pos{0.f, 0.f, 0.f};
    const scalar time{0.f};
    const vector3 mom{1.f * unit<scalar>::GeV, 0.f, 1.f * unit<scalar>::GeV};

    // Make a helix
    detray::detail::helix<traccc::default_algebra> hlx(
        pos, time, vector::normalize(mom), q / vector::norm(mom), B);

    // Make three spacepoints with the helix
    edm::measurement_collection::host measurements(host_mr);
    edm::spacepoint_collection::host spacepoints{host_mr};
    measurements.resize(3);
    spacepoints.reserve(3);
    spacepoints.push_back(
        {0, traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         traccc::utils::to_float_array<traccc::default_algebra>(
             hlx(50 * unit<scalar>::mm)),
         0.f, 0.f});
    spacepoints.push_back(
        {1, traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         traccc::utils::to_float_array<traccc::default_algebra>(
             hlx(100 * unit<scalar>::mm)),
         0.f, 0.f});
    spacepoints.push_back(
        {2, traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         traccc::utils::to_float_array<traccc::default_algebra>(
             hlx(150 * unit<scalar>::mm)),
         0.f, 0.f});

    // Make a seed from the three spacepoints
    edm::seed_collection::host seeds{host_mr};
    seeds.push_back({0, 1, 2, 0.0f});

    // Run track parameter estimation
    traccc::track_params_estimation_config track_params_estimation_config;
    traccc::host::track_params_estimation tp(track_params_estimation_config,
                                             host_mr);
    auto bound_params =
        tp(vecmem::get_data(measurements), vecmem::get_data(spacepoints),
           vecmem::get_data(seeds), B);

    // Make sure that the reconstructed momentum is equal to the original
    // momentum
    ASSERT_EQ(bound_params.size(), 1u);
    ASSERT_NEAR(bound_params[0].p(q), vector::norm(mom), 2.f * 1e-4);
    ASSERT_TRUE(bound_params[0].qop() < 0.f);
}

TEST(track_params_estimation, helix_positive_charge) {

    // Set B field
    const vector3 B{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                    2.f * unit<scalar>::T};

    // Track property
    const scalar q{1.f * unit<scalar>::e};
    const point3 pos{0.f, 0.f, 0.f};
    const scalar time{0.f};
    const vector3 mom{1.f * unit<scalar>::GeV, 0.f, 1.f * unit<scalar>::GeV};

    // Make a helix
    detray::detail::helix<traccc::default_algebra> hlx(
        pos, time, vector::normalize(mom), q / vector::norm(mom), B);

    // Make three spacepoints with the helix
    edm::measurement_collection::host measurements(host_mr);
    edm::spacepoint_collection::host spacepoints{host_mr};
    measurements.resize(3);
    spacepoints.reserve(3);
    spacepoints.push_back(
        {0, traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         traccc::utils::to_float_array<traccc::default_algebra>(
             hlx(50 * unit<scalar>::mm)),
         0.f, 0.f});
    spacepoints.push_back(
        {1, traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         traccc::utils::to_float_array<traccc::default_algebra>(
             hlx(100 * unit<scalar>::mm)),
         0.f, 0.f});
    spacepoints.push_back(
        {2, traccc::edm::spacepoint_collection::host::INVALID_MEASUREMENT_INDEX,
         traccc::utils::to_float_array<traccc::default_algebra>(
             hlx(150 * unit<scalar>::mm)),
         0.f, 0.f});

    // Make a seed from the three spacepoints
    edm::seed_collection::host seeds{host_mr};
    seeds.push_back({0, 1, 2, 0.0f});

    // Run track parameter estimation
    traccc::track_params_estimation_config track_params_estimation_config;
    traccc::host::track_params_estimation tp(track_params_estimation_config,
                                             host_mr);
    auto bound_params =
        tp(vecmem::get_data(measurements), vecmem::get_data(spacepoints),
           vecmem::get_data(seeds), B);

    // Make sure that the reconstructed momentum is equal to the original
    // momentum
    ASSERT_EQ(bound_params.size(), 1u);
    ASSERT_NEAR(bound_params[0].p(q), vector::norm(mom), 2.f * 1e-4);
    ASSERT_TRUE(bound_params[0].qop() > 0.f);
}
