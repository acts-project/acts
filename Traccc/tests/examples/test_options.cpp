/**
 * TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/options/details/value_array.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <sstream>
#include <stdexcept>

TEST(options, value_array_int) {

    traccc::opts::value_array<int, 3> parsed_vars;

    std::istringstream I("100:200:300");
    I >> parsed_vars;

    EXPECT_EQ(parsed_vars[0], 100);
    EXPECT_EQ(parsed_vars[1], 200);
    EXPECT_EQ(parsed_vars[2], 300);

    EXPECT_THROW(std::istringstream("100:200") >> parsed_vars,
                 std::invalid_argument);
    EXPECT_THROW(std::istringstream("100:200:300:400") >> parsed_vars,
                 std::invalid_argument);
}

TEST(options, value_array_float) {

    traccc::opts::value_array<float, 3> parsed_vars;

    std::istringstream I("1.1:2.2:3.3");
    I >> parsed_vars;

    EXPECT_FLOAT_EQ(parsed_vars[0], 1.1f);
    EXPECT_FLOAT_EQ(parsed_vars[1], 2.2f);
    EXPECT_FLOAT_EQ(parsed_vars[2], 3.3f);

    EXPECT_THROW(std::istringstream("1.1:2.2") >> parsed_vars,
                 std::invalid_argument);
    EXPECT_THROW(std::istringstream("1.1:2.2:3.3:4.4") >> parsed_vars,
                 std::invalid_argument);
}

TEST(options, value_array_double) {

    traccc::opts::value_array<double, 3> parsed_vars;

    std::istringstream I("1.1:2.2:3.3");
    I >> parsed_vars;

    EXPECT_DOUBLE_EQ(parsed_vars[0], 1.1);
    EXPECT_DOUBLE_EQ(parsed_vars[1], 2.2);
    EXPECT_DOUBLE_EQ(parsed_vars[2], 3.3);

    EXPECT_THROW(std::istringstream("1.1:2.2") >> parsed_vars,
                 std::invalid_argument);
    EXPECT_THROW(std::istringstream("1.1:2.2:3.3:4.4") >> parsed_vars,
                 std::invalid_argument);
}
