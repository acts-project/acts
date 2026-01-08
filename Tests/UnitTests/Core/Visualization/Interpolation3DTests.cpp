// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Visualization/Interpolation3D.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <numbers>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(VisualizationSuite)

BOOST_AUTO_TEST_CASE(SplineInterpolationEigen) {
  /// Define the input vector
  double R = 10.;
  std::vector<Acts::Vector3> inputs;

  // Interpolate the points options
  std::vector<Acts::Vector3> trajectory;

  // Empty in empty out check
  trajectory = Acts::Interpolation3D::spline(inputs, 10);
  BOOST_CHECK(trajectory.empty());

  for (double phi = 0; phi < 2 * std::numbers::pi;
       phi += std::numbers::pi / 4) {
    inputs.push_back(Acts::Vector3(R * std::cos(phi), R * std::sin(phi), 0.));
  }

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
    // Verify points remain in the XY plane
    CHECK_CLOSE_ABS(point.z(), 0., 0.1);
  }
}

BOOST_AUTO_TEST_CASE(SplineInterpolationArray) {
  /// Define the input vector
  std::vector<std::array<double, 3u>> inputs;

  for (double x = 0; x < 10; x += 1) {
    inputs.push_back({x, x * x, 0.});
  }

  // This time we keep the original hits
  auto trajectory = Acts::Interpolation3D::spline(inputs, 100, true);

  // Check the output type is correct
  constexpr bool isOutput =
      std::is_same_v<decltype(trajectory), decltype(inputs)>;
  BOOST_CHECK(isOutput);

  // Check the output size is correct
  BOOST_CHECK_EQUAL(trajectory.size(), 108);
}

BOOST_AUTO_TEST_CASE(SplineInterpolationErrors) {
  std::vector<std::array<double, 3u>> inputs;

  // Test with single point
  inputs.push_back({0., 0., 0.});
  auto result = Acts::Interpolation3D::spline(inputs, 10);
  BOOST_CHECK_EQUAL(result.size(), 1);

  // Test with two points
  inputs.push_back({1., 1., 1.});
  result = Acts::Interpolation3D::spline(inputs, 10);
  BOOST_CHECK_EQUAL(result.size(), 2);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
