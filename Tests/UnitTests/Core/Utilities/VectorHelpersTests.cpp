// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <limits>
#include <numbers>

using namespace Acts;
using namespace Acts::VectorHelpers;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(PhiFromVector) {
  // Test phi calculation for Eigen vectors
  CHECK_CLOSE_ABS(0.0, phi(Vector2{1, 0}), 1e-6);
  CHECK_CLOSE_ABS(std::numbers::pi / 2, phi(Vector2{0, 1}), 1e-6);
  CHECK_CLOSE_ABS(std::numbers::pi, phi(Vector2{-1, 0}), 1e-6);
  CHECK_CLOSE_ABS(-std::numbers::pi / 2, phi(Vector2{0, -1}), 1e-6);
  CHECK_CLOSE_ABS(std::numbers::pi / 4, phi(Vector2{1, 1}), 1e-6);

  // Test with 3D vectors
  CHECK_CLOSE_ABS(0.0, phi(Vector3{1, 0, 5}), 1e-6);
  CHECK_CLOSE_ABS(std::numbers::pi / 2, phi(Vector3{0, 1, -2}), 1e-6);

  // Test with dynamic vectors
  DynamicVector dynVec2(2);
  dynVec2 << 1, 1;
  CHECK_CLOSE_ABS(std::numbers::pi / 4, phi(dynVec2), 1e-6);

  DynamicVector dynVec3(3);
  dynVec3 << 0, 1, 5;
  CHECK_CLOSE_ABS(std::numbers::pi / 2, phi(dynVec3), 1e-6);
}

BOOST_AUTO_TEST_CASE(PhiFromObject) {
  // Test phi calculation for objects with phi() method
  struct MockObject {
    double phi() const { return 1.5; }
  };

  MockObject obj;
  CHECK_CLOSE_ABS(1.5, phi(obj), 1e-6);
}

BOOST_AUTO_TEST_CASE(PerpFromVector) {
  // Test transverse radius calculation
  CHECK_CLOSE_ABS(1.0, perp(Vector2{1, 0}), 1e-6);
  CHECK_CLOSE_ABS(1.0, perp(Vector2{0, 1}), 1e-6);
  CHECK_CLOSE_ABS(std::sqrt(2.0), perp(Vector2{1, 1}), 1e-6);
  CHECK_CLOSE_ABS(5.0, perp(Vector2{3, 4}), 1e-6);

  // Test with 3D vectors (should ignore z component)
  CHECK_CLOSE_ABS(1.0, perp(Vector3{1, 0, 10}), 1e-6);
  CHECK_CLOSE_ABS(5.0, perp(Vector3{3, 4, 100}), 1e-6);

  // Test with dynamic vectors
  DynamicVector dynVec(3);
  dynVec << 3, 4, 5;
  CHECK_CLOSE_ABS(5.0, perp(dynVec), 1e-6);
}

BOOST_AUTO_TEST_CASE(ThetaFromVector) {
  // Test theta calculation
  CHECK_CLOSE_ABS(std::numbers::pi / 2, theta(Vector3{1, 0, 0}), 1e-6);
  CHECK_CLOSE_ABS(std::numbers::pi / 2, theta(Vector3{0, 1, 0}), 1e-6);
  CHECK_CLOSE_ABS(0.0, theta(Vector3{0, 0, 1}), 1e-6);
  CHECK_CLOSE_ABS(std::numbers::pi, theta(Vector3{0, 0, -1}), 1e-6);
  CHECK_CLOSE_ABS(std::numbers::pi / 4, theta(Vector3{1, 0, 1}), 1e-6);
  CHECK_CLOSE_ABS(3 * std::numbers::pi / 4, theta(Vector3{1, 0, -1}), 1e-6);

  // Test with dynamic vectors
  DynamicVector dynVec(3);
  dynVec << 1, 0, 1;
  CHECK_CLOSE_ABS(std::numbers::pi / 4, theta(dynVec), 1e-6);

  dynVec << 0, 0, 1;
  CHECK_CLOSE_ABS(0.0, theta(dynVec), 1e-6);
}

BOOST_AUTO_TEST_CASE(EtaFromVector) {
  // Test pseudorapidity calculation
  CHECK_CLOSE_ABS(0.0, eta(Vector3{1, 0, 0}), 1e-6);
  CHECK_CLOSE_ABS(0.0, eta(Vector3{0, 1, 0}), 1e-6);

  // Test special case: vector parallel to z-axis
  BOOST_CHECK_EQUAL(eta(Vector3{0, 0, 1}),
                    std::numeric_limits<double>::infinity());
  BOOST_CHECK_EQUAL(eta(Vector3{0, 0, -1}),
                    -std::numeric_limits<double>::infinity());

  // Test finite values
  CHECK_CLOSE_ABS(std::asinh(1.0), eta(Vector3{1, 0, 1}), 1e-6);
  CHECK_CLOSE_ABS(std::asinh(-1.0), eta(Vector3{1, 0, -1}), 1e-6);
  CHECK_CLOSE_ABS(std::asinh(2.0), eta(Vector3{1, 1, 2 * std::sqrt(2)}), 1e-6);

  // Test with dynamic vectors
  DynamicVector dynVec(3);
  dynVec << 1, 0, 1;
  CHECK_CLOSE_ABS(std::asinh(1.0), eta(dynVec), 1e-6);

  dynVec << 0, 0, 1;
  BOOST_CHECK_EQUAL(eta(dynVec), std::numeric_limits<double>::infinity());

  dynVec << 0, 0, -1;
  BOOST_CHECK_EQUAL(eta(dynVec), -std::numeric_limits<double>::infinity());
}

BOOST_AUTO_TEST_CASE(EvaluateTrigonomics) {
  // Test trigonometric evaluation
  Vector3 dir{1, 0, 0};  // pointing in x direction
  auto trig = evaluateTrigonomics(dir);

  CHECK_CLOSE_ABS(1.0, trig[0], 1e-6);  // cos(phi)
  CHECK_CLOSE_ABS(0.0, trig[1], 1e-6);  // sin(phi)
  CHECK_CLOSE_ABS(0.0, trig[2], 1e-6);  // cos(theta)
  CHECK_CLOSE_ABS(1.0, trig[3], 1e-6);  // sin(theta)

  // Test with y direction
  dir = Vector3{0, 1, 0};
  trig = evaluateTrigonomics(dir);

  CHECK_CLOSE_ABS(0.0, trig[0], 1e-6);  // cos(phi)
  CHECK_CLOSE_ABS(1.0, trig[1], 1e-6);  // sin(phi)
  CHECK_CLOSE_ABS(0.0, trig[2], 1e-6);  // cos(theta)
  CHECK_CLOSE_ABS(1.0, trig[3], 1e-6);  // sin(theta)

  // Test with diagonal direction (avoid z-axis singularity)
  dir = Vector3{1, 1, 1};
  dir.normalize();
  trig = evaluateTrigonomics(dir);

  // For normalized {1,1,1}, cos(theta) = 1/sqrt(3), sin(theta) = sqrt(2/3)
  // cos(phi) = sin(phi) = 1/sqrt(2) (45 degree angle in xy plane)
  CHECK_CLOSE_ABS(1.0 / std::sqrt(2), trig[0], 1e-6);    // cos(phi)
  CHECK_CLOSE_ABS(1.0 / std::sqrt(2), trig[1], 1e-6);    // sin(phi)
  CHECK_CLOSE_ABS(1.0 / std::sqrt(3), trig[2], 1e-6);    // cos(theta)
  CHECK_CLOSE_ABS(std::sqrt(2.0 / 3.0), trig[3], 1e-6);  // sin(theta)
}

BOOST_AUTO_TEST_CASE(CastAxisDirection) {
  using enum AxisDirection;
  Vector3 pos{3, 4, 5};

  // Test all axis directions
  CHECK_CLOSE_ABS(3.0, cast(pos, AxisX), 1e-6);
  CHECK_CLOSE_ABS(4.0, cast(pos, AxisY), 1e-6);
  CHECK_CLOSE_ABS(5.0, cast(pos, AxisZ), 1e-6);
  CHECK_CLOSE_ABS(5.0, cast(pos, AxisR), 1e-6);  // sqrt(3²+4²)
  CHECK_CLOSE_ABS(std::atan2(4, 3), cast(pos, AxisPhi), 1e-6);
  CHECK_CLOSE_ABS(5.0 * std::atan2(4, 3), cast(pos, AxisRPhi), 1e-6);
  CHECK_CLOSE_ABS(std::atan2(5, 5), cast(pos, AxisTheta), 1e-6);
  CHECK_CLOSE_ABS(std::asinh(1.0), cast(pos, AxisEta), 1e-6);
  CHECK_CLOSE_ABS(std::sqrt(50), cast(pos, AxisMag), 1e-6);
}

BOOST_AUTO_TEST_CASE(CrossProduct) {
  SquareMatrix3 m;
  m << 1, 0, 0, 0, 1, 0, 0, 0, 1;

  Vector3 v{1, 0, 0};

  SquareMatrix3 result = cross(m, v);

  // Check that each column is the cross product of the corresponding column of
  // m with v
  Vector3 expected_col0 = m.col(0).cross(v);  // [1,0,0] x [1,0,0] = [0,0,0]
  Vector3 expected_col1 = m.col(1).cross(v);  // [0,1,0] x [1,0,0] = [0,0,-1]
  Vector3 expected_col2 = m.col(2).cross(v);  // [0,0,1] x [1,0,0] = [0,1,0]

  CHECK_CLOSE_ABS(expected_col0[0], result.col(0)[0], 1e-6);
  CHECK_CLOSE_ABS(expected_col0[1], result.col(0)[1], 1e-6);
  CHECK_CLOSE_ABS(expected_col0[2], result.col(0)[2], 1e-6);

  CHECK_CLOSE_ABS(expected_col1[0], result.col(1)[0], 1e-6);
  CHECK_CLOSE_ABS(expected_col1[1], result.col(1)[1], 1e-6);
  CHECK_CLOSE_ABS(expected_col1[2], result.col(1)[2], 1e-6);

  CHECK_CLOSE_ABS(expected_col2[0], result.col(2)[0], 1e-6);
  CHECK_CLOSE_ABS(expected_col2[1], result.col(2)[1], 1e-6);
  CHECK_CLOSE_ABS(expected_col2[2], result.col(2)[2], 1e-6);
}

BOOST_AUTO_TEST_CASE(PositionFromVector4) {
  Vector4 pos4{1, 2, 3, 4};
  auto pos3 = position(pos4);

  CHECK_CLOSE_ABS(1.0, pos3[0], 1e-6);
  CHECK_CLOSE_ABS(2.0, pos3[1], 1e-6);
  CHECK_CLOSE_ABS(3.0, pos3[2], 1e-6);
}

BOOST_AUTO_TEST_CASE(PositionFromFreeVector) {
  FreeVector params = FreeVector::Zero();
  params[eFreePos0] = 1;
  params[eFreePos1] = 2;
  params[eFreePos2] = 3;

  auto pos3 = position(params);

  CHECK_CLOSE_ABS(1.0, pos3[0], 1e-6);
  CHECK_CLOSE_ABS(2.0, pos3[1], 1e-6);
  CHECK_CLOSE_ABS(3.0, pos3[2], 1e-6);
}

BOOST_AUTO_TEST_CASE(MakeVector4) {
  Vector3 vec3{1, 2, 3};
  double w = 4.0;

  auto vec4 = makeVector4(vec3, w);

  CHECK_CLOSE_ABS(1.0, vec4[ePos0], 1e-6);
  CHECK_CLOSE_ABS(2.0, vec4[ePos1], 1e-6);
  CHECK_CLOSE_ABS(3.0, vec4[ePos2], 1e-6);
  CHECK_CLOSE_ABS(4.0, vec4[eTime], 1e-6);
}

BOOST_AUTO_TEST_CASE(IncidentAngles) {
  Vector3 direction{1, 1, 1};
  direction.normalize();

  // Identity rotation (global == local)
  RotationMatrix3 globalToLocal = RotationMatrix3::Identity();

  auto angles = incidentAngles(direction, globalToLocal);

  // In local frame, direction is still {1,1,1}/sqrt(3)
  // phi = atan2(z, x) = atan2(1/sqrt(3), 1/sqrt(3)) = pi/4
  // theta = atan2(z, y) = atan2(1/sqrt(3), 1/sqrt(3)) = pi/4
  CHECK_CLOSE_ABS(std::numbers::pi / 4, angles.first, 1e-6);
  CHECK_CLOSE_ABS(std::numbers::pi / 4, angles.second, 1e-6);

  // Test with rotation around z-axis by 90 degrees
  RotationMatrix3 rotZ;
  rotZ << 0, 1, 0, -1, 0, 0, 0, 0, 1;

  angles = incidentAngles(Vector3{1, 0, 0}, rotZ);
  // After rotation: {1,0,0} -> {0,-1,0}
  // phi = atan2(0, 0) = 0, theta = atan2(0, -1) = pi (atan2 of 0,-1)
  CHECK_CLOSE_ABS(0.0, angles.first, 1e-6);
  CHECK_CLOSE_ABS(std::numbers::pi, angles.second, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
