/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/geometry/module_map.hpp"
#include "traccc/io/utils.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <map>
#include <string>

/*
 * Simple test of this map using integers and strings.
 */
TEST(geometry, module_map_simple) {
    std::map<std::size_t, std::string> inp{{0, "zero"},    {1, "one"},
                                           {2, "two"},     {5, "five"},
                                           {11, "eleven"}, {21, "twenty-one"}};

    /*
     * Convert the old-style map to a module map.
     */
    traccc::module_map<std::size_t, std::string> map(inp);

    /*
     * Check whether the two maps are the same size!
     */
    ASSERT_EQ(map.size(), inp.size());

    ASSERT_EQ(map.at(0), "zero");
    ASSERT_EQ(map.at(1), "one");
    ASSERT_EQ(map.at(2), "two");
    ASSERT_EQ(map.at(5), "five");
    ASSERT_EQ(map.at(11), "eleven");
    ASSERT_EQ(map.at(21), "twenty-one");
}

/*
 * Ensure at and operator[] do the same thing.
 */
TEST(geometry, module_map_operator_eq) {
    std::map<std::size_t, std::string> inp{{0, "zero"},    {1, "one"},
                                           {2, "two"},     {5, "five"},
                                           {11, "eleven"}, {21, "twenty-one"}};

    /*
     * Convert the old-style map to a module map.
     */
    traccc::module_map<std::size_t, std::string> map(inp);

    /*
     * Check if these are all the same.
     */
    ASSERT_EQ(map.at(0), map[0]);
    ASSERT_EQ(map.at(1), map[1]);
    ASSERT_EQ(map.at(2), map[2]);
    ASSERT_EQ(map.at(5), map[5]);
    ASSERT_EQ(map.at(11), map[11]);
    ASSERT_EQ(map.at(21), map[21]);
}

/*
 * Sanity check for the empty method.
 */
TEST(geometry, module_map_empty) {
    std::map<std::size_t, std::string> inp{{0, "zero"},    {1, "one"},
                                           {2, "two"},     {5, "five"},
                                           {11, "eleven"}, {21, "twenty-one"}};

    /*
     * Convert the old-style map to a module map.
     */
    traccc::module_map<std::size_t, std::string> map(inp);

    ASSERT_FALSE(map.empty());
}

/*
 * Check whether the contains method works.
 */
TEST(geometry, module_map_contains) {
    std::map<std::size_t, std::string> inp{{0, "zero"},    {1, "one"},
                                           {2, "two"},     {5, "five"},
                                           {11, "eleven"}, {21, "twenty-one"}};

    /*
     * Convert the old-style map to a module map.
     */
    traccc::module_map<std::size_t, std::string> map(inp);

    ASSERT_TRUE(map.contains(0));
    ASSERT_TRUE(map.contains(1));
    ASSERT_TRUE(map.contains(2));
    ASSERT_TRUE(map.contains(5));
    ASSERT_TRUE(map.contains(11));
    ASSERT_TRUE(map.contains(21));
    ASSERT_FALSE(map.contains(3));
    ASSERT_FALSE(map.contains(4));
    ASSERT_FALSE(map.contains(6));
    ASSERT_FALSE(map.contains(15));
    ASSERT_FALSE(map.contains(510482));
}

/*
 * Check if the map fails correctly.
 */
TEST(geometry, module_map_failure) {
    std::map<std::size_t, std::string> inp{{0, "zero"}, {1, "one"}, {2, "two"}};

    /*
     * Convert the old-style map to a module map.
     */
    traccc::module_map<std::size_t, std::string> map(inp);

    /*
     * This value is not in the map, so it should throw.
     */
    ASSERT_THROW(map.at(100), std::out_of_range);
}
