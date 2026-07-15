/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "traccc/definitions/primitives.hpp"
#include "traccc/seeding/grids/axis.hpp"
#include "traccc/seeding/grids/grid2.hpp"
#include "traccc/seeding/grids/populator.hpp"
#include "traccc/seeding/grids/serializer2.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace traccc;

GTEST_TEST(traccc_grid2, grid2_replace_populator) {
    vecmem::host_memory_resource host_mr;

    serializer2 serializer;

    using grid2r = grid2<replace_populator, axis2::regular, axis2::regular,
                         decltype(serializer)>;
    typename grid2r::axis_p0_type xaxis{10u, -5.f, 5.f, host_mr};
    typename grid2r::axis_p1_type yaxis{10u, -5.f, 5.f, host_mr};

    grid2r g2(std::move(xaxis), std::move(yaxis), host_mr);

    // Test the initialization
    point2 p = {-4.5f, -4.5f};
    for (unsigned int ib0 = 0u; ib0 < 10u; ++ib0) {
        for (unsigned int ib1 = 0u; ib1 < 10u; ++ib1) {
            p = {-4.5f + static_cast<scalar>(ib0),
                 -4.5f + static_cast<scalar>(ib1)};
            EXPECT_EQ(g2.bin(p), std::numeric_limits<unsigned int>::max());
        }
    }

    p = {-4.5f, -4.5f};
    // Fill and read
    g2.populate(p, 3u);
    EXPECT_EQ(g2.bin(p), 3u);

    // Fill and read two times, fill first 0-99, then 100-199
    for (unsigned int il = 0u; il < 2u; ++il) {
        unsigned int counter = il * 100u;
        for (unsigned int ib0 = 0u; ib0 < 10u; ++ib0) {
            for (unsigned int ib1 = 0; ib1 < 10u; ++ib1) {
                p = {-4.5f + static_cast<scalar>(ib0),
                     -4.5f + static_cast<scalar>(ib1)};
                g2.populate(p, counter + 0u);
                EXPECT_EQ(g2.bin(p), counter++);
            }
        }
    }

    // A zone test w/o neighbour hood
    p = {-4.5f, -4.5f};
    auto test = g2.zone(p);
    vecmem::vector<unsigned int> expect = {100u};
    EXPECT_EQ(test, expect);

    // A zone test with neighbour hood
    p = {0.5f, 0.5f};

    std::array<unsigned int, 2> zone11 = {1u, 1u};
    std::array<unsigned int, 2> zone22 = {2u, 2u};

    test = g2.zone(p, {zone11, zone22}, true);
    expect = {143u, 144u, 145u, 146u, 147u, 153u, 154u, 155u,
              156u, 157u, 163u, 164u, 165u, 166u, 167u};
    EXPECT_EQ(test, expect);

    using grid2cc = grid2<replace_populator, axis2::circular, axis2::regular,
                          decltype(serializer)>;

    typename grid2cc::axis_p0_type circular{4u, -2.f, 2.f, host_mr};
    typename grid2cc::axis_p1_type closed{5u, 0.f, 5.f, host_mr};

    grid2cc g2cc(std::move(circular), std::move(closed), host_mr);
    unsigned int counter = 0u;
    for (unsigned icl = 0u; icl < 5u; ++icl) {
        for (unsigned ici = 0u; ici < 4u; ++ici) {
            p = {-1.5f + static_cast<scalar>(ici),
                 0.5f + static_cast<scalar>(icl)};
            g2cc.populate(p, counter++);
        }
    }

    // A zone test for circular testing
    p = {1.5f, 2.5f};
    test = g2cc.zone(p, {zone11, zone11}, true);
    expect = {4u, 6u, 7u, 8u, 10u, 11u, 12u, 14u, 15u};
    EXPECT_EQ(test, expect);
}

GTEST_TEST(traccc_grid2, grid2_complete_populator) {
    constexpr auto invalid_index{std::numeric_limits<unsigned int>::max()};

    vecmem::host_memory_resource host_mr;

    serializer2 serializer;

    using grid2r =
        grid2<complete_populator, axis2::regular, axis2::regular,
              decltype(serializer), vecmem::vector, vecmem::jagged_vector,
              std::array, detray::tuple, unsigned int, false, 3>;

    typename grid2r::axis_p0_type xaxis{2u, -1.f, 1.f, host_mr};
    typename grid2r::axis_p1_type yaxis{2u, -1.f, 1.f, host_mr};

    grid2r g2(std::move(xaxis), std::move(yaxis), host_mr);

    // Test the initialization
    point2 p = {-0.5f, -0.5f};
    grid2r::populator_type::store_value invalid = {invalid_index, invalid_index,
                                                   invalid_index};
    for (unsigned int ib0 = 0u; ib0 < 2u; ++ib0) {
        for (unsigned int ib1 = 0u; ib1 < 2u; ++ib1) {
            p = {-0.5f + static_cast<scalar>(ib0),
                 -0.5f + static_cast<scalar>(ib1)};
            EXPECT_EQ(g2.bin(p), invalid);
        }
    }

    // Fill and read
    p = {-0.5f, -0.5f};
    g2.populate(p, 4u);

    grid2r::populator_type::store_value expected = {4u, invalid_index,
                                                    invalid_index};
    auto test = g2.bin(p);
    EXPECT_EQ(test, expected);

    auto zone_test = g2.zone(p);
    vecmem::vector<unsigned int> zone_expected = {4u};
    EXPECT_EQ(zone_test, zone_expected);

    g2.populate(p, 2u);
    expected = {4u, 2u, invalid_index};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    g2.populate(p, 7u);
    expected = {4u, 2u, 7u};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    // Bin is completed, new entry is ignored
    g2.populate(p, 16u);
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    std::array<unsigned int, 2> zone00 = {0u, 0u};
    std::array<unsigned int, 2> zone11 = {1u, 1u};

    // Zone test of a complete bin
    zone_test = g2.zone(p, {zone00, zone00});
    zone_expected = {4u, 2u, 7u};
    EXPECT_EQ(zone_test, zone_expected);

    // Fill some other bins
    p = {0.5f, -0.5f};
    g2.populate(p, 16u);

    p = {0.5f, 0.5f};
    g2.populate(p, 17u);
    g2.populate(p, 18u);

    zone_test = g2.zone(p, {zone11, zone11});
    zone_expected = {4u, 2u, 7u, 16u, 17u, 18u};
    EXPECT_EQ(zone_test, zone_expected);
}

GTEST_TEST(traccc_grid2, grid2_attach_populator) {
    vecmem::host_memory_resource host_mr;

    serializer2 serializer;

    using grid2r = grid2<attach_populator, axis2::regular, axis2::regular,
                         decltype(serializer)>;
    typename grid2r::axis_p0_type xaxis{2u, -1.f, 1.f, host_mr};
    typename grid2r::axis_p1_type yaxis{2u, -1.f, 1.f, host_mr};

    grid2r g2(std::move(xaxis), std::move(yaxis), host_mr);

    // Test the initialization
    point2 p = {-0.5f, -0.5f};
    grid2r::populator_type::store_value invalid = {};
    for (unsigned int ib0 = 0u; ib0 < 2u; ++ib0) {
        for (unsigned int ib1 = 0u; ib1 < 2u; ++ib1) {
            p = {-0.5f + static_cast<scalar>(ib0),
                 -0.5f + static_cast<scalar>(ib1)};
            EXPECT_EQ(g2.bin(p), invalid);
        }
    }

    p = {-0.5f, -0.5f};
    g2.populate(p, 4u);

    grid2r::populator_type::store_value expected = {4u};
    auto test = g2.bin(p);
    EXPECT_EQ(test, expected);

    auto zone_test = g2.zone(p);
    vecmem::vector<unsigned int> zone_expected = {4u};
    EXPECT_EQ(zone_test, zone_expected);

    p = {-0.5f, 0.5f};
    g2.populate(p, 9u);

    p = {0.5f, -0.5f};
    g2.populate(p, 1u);

    p = {0.5f, 0.5f};
    g2.populate(p, 7u);

    expected = {7u};
    test = g2.bin(p);
    EXPECT_EQ(test, expected);

    std::array<unsigned int, 2> zone11 = {1u, 1u};

    zone_test = g2.zone(p, {zone11, zone11}, true);
    zone_expected = {1u, 4u, 7u, 9u};
    EXPECT_EQ(zone_test, zone_expected);
}

GTEST_TEST(traccc_grid2, grid2_shift) {
    vecmem::host_memory_resource host_mr;

    serializer2 serializer;

    using grid2r = grid2<replace_populator, axis2::regular, axis2::regular,
                         decltype(serializer)>;

    typename grid2r::axis_p0_type xaxis{10u, -5.f, 5.f, host_mr};
    typename grid2r::axis_p1_type yaxis{10u, -5.f, 5.f, host_mr};

    grid2r g2(std::move(xaxis), std::move(yaxis), host_mr, 0);

    // Test the initialization
    point2 p = {-4.5f, -4.5f};
    EXPECT_EQ(g2.bin(p), 0u);

    g2.shift(8u);
    EXPECT_EQ(g2.bin(p), 8u);
}

GTEST_TEST(traccc_grid2, grid2_irregular_replace) {
    vecmem::host_memory_resource host_mr;

    replace_populator<> replacer;
    serializer2 serializer;

    using grid2ir = grid2<replace_populator, axis2::irregular, axis2::irregular,
                          decltype(serializer)>;

    typename grid2ir::axis_p0_type xaxis{
        {-3.f, -2.f, 1.f, 0.5f, 0.7f, 0.71f, 4.f, 1000.f}, host_mr};
    typename grid2ir::axis_p1_type yaxis{{0.1f, 0.8f, 0.9f, 10.f, 12.f, 15.f},
                                         host_mr};

    grid2ir g2(std::move(xaxis), std::move(yaxis), host_mr);

    point2 p = {-0.5f, 0.5f};
    g2.populate(p, 4u);
    EXPECT_EQ(g2.bin(p), 4u);
}
