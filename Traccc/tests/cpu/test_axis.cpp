/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/seeding/grids/axis.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

// detray core

using namespace traccc;

GTEST_TEST(traccc_grid2, regular_closed_axis) {
    vecmem::host_memory_resource resource;

    axis2::regular<> ten_bins{10u, -3.f, 7.f, resource};
    // N bins
    EXPECT_EQ(ten_bins.bins(), 10u);
    // Axis bin access
    EXPECT_EQ(ten_bins.bin(-4.f), 0u);
    EXPECT_EQ(ten_bins.bin(2.5f), 5u);
    EXPECT_EQ(ten_bins.bin(8.f), 9u);

    // Axis range access - binned (symmetric & asymmetric)
    std::array<unsigned int, 2> zone00 = {0u, 0u};
    std::array<unsigned int, 2> zone01 = {0u, 1u};
    std::array<unsigned int, 2> zone11 = {1u, 1u};
    std::array<unsigned int, 2> zone44 = {4u, 4u};
    std::array<unsigned int, 2> zone55 = {5u, 5u};

    std::array<unsigned int, 2> expected_range = {5u, 5u};
    EXPECT_EQ(ten_bins.range(2.5f, zone00), expected_range);
    expected_range = {4u, 6u};
    EXPECT_EQ(ten_bins.range(2.5f, zone11), expected_range);
    expected_range = {5u, 6u};
    EXPECT_EQ(ten_bins.range(2.5f, zone01), expected_range);
    expected_range = {0u, 8u};
    EXPECT_EQ(ten_bins.range(1.5f, zone44), expected_range);
    expected_range = {3u, 9u};
    EXPECT_EQ(ten_bins.range(5.5f, zone55), expected_range);

    // Axis sequence access - binned (symmetric & asymmetric)
    vecmem::vector<unsigned int> expected_zone = {5u};
    EXPECT_EQ(ten_bins.zone(2.5f, zone00), expected_zone);
    expected_zone = {5u, 6u};
    EXPECT_EQ(ten_bins.zone(2.5f, zone01), expected_zone);
    expected_zone = {4u, 5u, 6u};
    EXPECT_EQ(ten_bins.zone(2.5f, zone11), expected_zone);
    expected_zone = {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u};
    EXPECT_EQ(ten_bins.zone(1.5f, zone44), expected_zone);

    // Axis range access - scalar (symmteric & asymmetric)
    std::array<scalar, 2> szone00 = {0.f, 0.f};
    std::array<scalar, 2> sepsilon = {0.01f, 0.01f};
    std::array<scalar, 2> szone11 = {1.f, 1.f};
    std::array<scalar, 2> szoneAll = {10.f, 10.f};

    expected_range = {5u, 5u};
    EXPECT_EQ(ten_bins.range(2.5f, szone00), expected_range);
    EXPECT_EQ(ten_bins.range(2.5f, sepsilon), expected_range);
    expected_range = {4u, 6u};
    EXPECT_EQ(ten_bins.range(2.5f, szone11), expected_range);
    expected_range = {0u, 9u};
    EXPECT_EQ(ten_bins.range(2.5f, szoneAll), expected_range);

    // Axis sequence acces - scalar (symmteric & asymmetric)
    expected_zone = {5u};
    EXPECT_EQ(ten_bins.zone(2.5f, szone00), expected_zone);
    EXPECT_EQ(ten_bins.zone(2.5f, sepsilon), expected_zone);
    expected_zone = {4u, 5u, 6u};
    EXPECT_EQ(ten_bins.zone(2.5f, szone11), expected_zone);
    expected_zone = {0u, 1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u};
    EXPECT_EQ(ten_bins.zone(2.5f, szoneAll), expected_zone);
}

GTEST_TEST(traccc_grid2, regular_circular_axis) {
    vecmem::host_memory_resource resource;

    constexpr scalar epsilon{10.f * std::numeric_limits<scalar>::epsilon()};

    // Let's say 36 modules, but with 4 directly at 0, pi/2, pi, -pi2
    scalar half_module{constant<scalar>::pi / 72.f};
    scalar phi_min = -constant<scalar>::pi + half_module;
    scalar phi_max = constant<scalar>::pi - half_module;
    axis2::circular<> full_pi = {36u, phi_min, phi_max, resource};
    // N bins
    EXPECT_EQ(full_pi.bins(), 36u);
    // Axis bin access
    EXPECT_EQ(full_pi.bin(constant<scalar>::pi - epsilon), 0u);
    EXPECT_EQ(full_pi.bin(constant<scalar>::pi + epsilon), 0u);
    EXPECT_EQ(full_pi.bin(0u), 18u);
    // Remap test
    EXPECT_EQ(full_pi.remap(4, -1), 3u);
    EXPECT_EQ(full_pi.remap(4, 1), 5u);
    EXPECT_EQ(full_pi.remap(0, -1), 35u);
    EXPECT_EQ(full_pi.remap(0, -2), 34u);
    EXPECT_EQ(full_pi.remap(1, -1), 0u);
    EXPECT_EQ(full_pi.remap(35, 1), 0u);

    // Axis range access - binned (symmetric & asymmetric)
    std::array<unsigned int, 2> zone00 = {0u, 0u};
    std::array<unsigned int, 2> zone01 = {0u, 1u};
    std::array<unsigned int, 2> zone11 = {1u, 1u};
    std::array<unsigned int, 2> zone22 = {2u, 2u};

    std::array<unsigned int, 2> expected_range = {0u, 0u};
    EXPECT_EQ(full_pi.range(constant<scalar>::pi + epsilon, zone00),
              expected_range);
    expected_range = {0u, 1u};
    EXPECT_EQ(full_pi.range(constant<scalar>::pi + epsilon, zone01),
              expected_range);
    expected_range = {35u, 1u};
    EXPECT_EQ(full_pi.range(constant<scalar>::pi + epsilon, zone11),
              expected_range);
    expected_range = {34u, 2u};
    EXPECT_EQ(full_pi.range(constant<scalar>::pi + epsilon, zone22),
              expected_range);

    // Zone test - binned
    vecmem::vector<unsigned int> expected_zone = {34u, 35u, 0u, 1u, 2u};
    EXPECT_EQ(full_pi.zone(constant<scalar>::pi + epsilon, zone22),
              expected_zone);

    // Axis range access - scalar (symmetric & asymmteric)
    std::array<scalar, 2> szone00 = {0.f, 0.f};
    std::array<scalar, 2> szoneEpsilon = {0.5f * epsilon, 0.5f * epsilon};
    scalar bin_step =
        (full_pi.max - full_pi.min) / static_cast<scalar>(full_pi.bins());
    std::array<scalar, 2> szone22 = {2.f * bin_step, 2.f * bin_step};

    expected_range = {0u, 0u};
    EXPECT_EQ(full_pi.range(constant<scalar>::pi + epsilon, szone00),
              expected_range);
    EXPECT_EQ(full_pi.range(constant<scalar>::pi + epsilon, szoneEpsilon),
              expected_range);

    expected_range = {34u, 2u};
    EXPECT_EQ(full_pi.range(constant<scalar>::pi + epsilon, szone22),
              expected_range);

    expected_zone = {34u, 35u, 0u, 1u, 2u};
    EXPECT_EQ(full_pi.zone(constant<scalar>::pi + epsilon, szone22),
              expected_zone);
}

GTEST_TEST(traccc_grid2, irregular_closed_axis) {
    vecmem::host_memory_resource resource;

    axis2::irregular<> nonreg({-3.f, 1.f, 2.f, 4.f, 8.f, 12.f}, resource);

    // Axis bin access
    //
    // N bins
    EXPECT_EQ(nonreg.bins(), 5u);
    // Bin tests
    EXPECT_EQ(nonreg.bin(-2), 0u);
    EXPECT_EQ(nonreg.bin(10), 4u);
    // Underflow test
    EXPECT_EQ(nonreg.bin(-4), 0u);
    // Overflow test
    EXPECT_EQ(nonreg.bin(14), 4u);

    // Axis range access - binned  (symmetric & asymmetric)
    std::array<unsigned int, 2> zone01 = {0u, 1u};
    std::array<unsigned int, 2> zone11 = {1u, 1u};
    std::array<unsigned int, 2> zone22 = {2u, 2u};

    std::array<unsigned int, 2> expected_range = {1u, 3u};
    EXPECT_EQ(nonreg.range(3.f, zone11), expected_range);

    expected_range = {2u, 3u};
    EXPECT_EQ(nonreg.range(3.f, zone01), expected_range);

    std::array<unsigned int, 2> expected_range_truncated_low = {0u, 1u};
    EXPECT_EQ(nonreg.range(0.f, zone11), expected_range_truncated_low);

    std::array<unsigned int, 2> expected_range_truncated_high = {2u, 4u};
    EXPECT_EQ(nonreg.range(10.f, zone22), expected_range_truncated_high);

    // Axis sequence access - binned
    vecmem::vector<unsigned int> expected_zone = {1u, 2u, 3u};
    EXPECT_EQ(nonreg.zone(3.f, zone11), expected_zone);

    vecmem::vector<unsigned int> expected_zone_truncated_low = {0u, 1u};
    EXPECT_EQ(nonreg.zone(0.f, zone11), expected_zone_truncated_low);

    vecmem::vector<unsigned int> expected_zone_truncated_high = {2u, 3u, 4u};
    EXPECT_EQ(nonreg.zone(10.f, zone22), expected_zone_truncated_high);

    // Axis range access - scalar
    std::array<scalar, 2> szone00 = {0.f, 0.f};
    std::array<scalar, 2> szone10 = {1.5f, 0.2f};
    expected_range = {2u, 2u};
    EXPECT_EQ(nonreg.range(3.f, szone00), expected_range);

    expected_range = {1u, 2u};
    EXPECT_EQ(nonreg.range(3.f, szone10), expected_range);

    // Axis sequence access - scalar
    expected_zone = {2u};
    EXPECT_EQ(nonreg.zone(3.f, szone00), expected_zone);

    expected_zone = {1u, 2u};
    EXPECT_EQ(nonreg.zone(3.f, szone10), expected_zone);
}
