/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/utils/ranges.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace traccc;

// Test eta_to_theta_range function
TEST(ranges, eta_to_theta) {

    // Comparing with the values provided by
    // https://www.star.bnl.gov/~dmitry/calc2.html
    ASSERT_NEAR(eta_to_theta(2.3f), 11.451f * traccc::unit<scalar>::degree,
                1e-3f * traccc::unit<scalar>::degree);
    ASSERT_NEAR(eta_to_theta(-0.87f), 134.537f * traccc::unit<scalar>::degree,
                1e-3f * traccc::unit<scalar>::degree);

    std::array<scalar, 2> eta_range{-2.44f, 3.13f};
    const auto theta_range = eta_to_theta_range(eta_range);

    ASSERT_NEAR(theta_range[0], 5.007f * traccc::unit<scalar>::degree,
                1e-3f * traccc::unit<scalar>::degree);
    ASSERT_NEAR(theta_range[1], 170.037f * traccc::unit<scalar>::degree,
                1e-3f * traccc::unit<scalar>::degree);
}

// Test theta_to_eta_range function
TEST(ranges, theta_to_eta) {

    // Comparing with the values provided by
    // https://www.star.bnl.gov/~dmitry/calc2.html
    ASSERT_NEAR(theta_to_eta(132.8f * traccc::unit<scalar>::degree), -0.828f,
                1e-3f);
    ASSERT_NEAR(theta_to_eta(90.f * traccc::unit<scalar>::degree), 0.f, 1e-3f);

    std::array<scalar, 2> theta_range{45.f * traccc::unit<scalar>::degree,
                                      175.f * traccc::unit<scalar>::degree};
    const auto eta_range = theta_to_eta_range(theta_range);

    ASSERT_NEAR(eta_range[0], -3.131f, 1e-3f);
    ASSERT_NEAR(eta_range[1], 0.881f, 1e-3f);
}
