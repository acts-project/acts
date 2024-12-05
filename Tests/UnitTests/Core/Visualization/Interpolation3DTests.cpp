// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <numbers>

#include "Acts/Visualization/Interpolation3D.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Visualization)

BOOST_AUTO_TEST_CASE(SplineInterpolation) {

    /// Define the input vector
    double R = 10.;
    std::vector<Acts::Vector3> inputs;
    for (double phi = 0; phi < 2 * std::numbers::pi; phi += std::numbers::pi / 4) {
        inputs.push_back(Acts::Vector3(R * cos(phi), R * sin(phi), 0.));
    }

    // Interpolate the points options
    std::vector<Acts::Vector3> trajectory;

    // (0) - nothing happens
    trajectory = Acts::Interpolation3D::spline(inputs, 1);
    // Check input and output size are the same
    BOOST_CHECK_EQUAL(trajectory.size(), inputs.size());

    // (1) - interpolate between the points with 12 points in total
    trajectory = Acts::Interpolation3D::spline(inputs, 12);
    // Check the output size is correct
    BOOST_CHECK_EQUAL(trajectory.size(), 12);

    for (const auto& point : trajectory) {
        // Check the interpolated points are on the circle
        // with a tolerance of course
        CHECK_CLOSE_ABS(point.norm(), R, 0.1);
    }

}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test