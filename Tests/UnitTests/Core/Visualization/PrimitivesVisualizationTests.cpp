// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Visualization/EventDataVisualization.hpp"
#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/IVisualization.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"
#include "Acts/Visualization/PlyVisualization.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(Visualization)

// Test on a plane
auto identity = std::make_shared<Transform3D>(Transform3D::Identity());
auto rectangle = std::make_shared<RectangleBounds>(10., 10.);
auto plane = Surface::makeShared<PlaneSurface>(identity, rectangle);

// Test context
GeometryContext gctx = GeometryContext();

/// The tests in this section are regression tests only in order
/// to catch any unexpected changes in the output format.
///
BOOST_AUTO_TEST_CASE(VisualizationHelpers) {
  // No correlation, fully summetric
  ActsSymMatrixD<2> covariance;
  covariance << 4., 0., 0., 4.;
  auto decops = Acts::Visualization::decomposeCovariance(covariance);
  BOOST_CHECK(decops[0] == 4.);
  BOOST_CHECK(decops[1] == 4.);
  BOOST_CHECK(decops[2] == 0.);

  // Fully positively correlated
  covariance.setZero();
  covariance << 4., 4., 4., 4.;
  decops = Acts::Visualization::decomposeCovariance(covariance);
  BOOST_CHECK(decops[0] == 8.);
  BOOST_CHECK(decops[1] == 0.);
  CHECK_CLOSE_ABS(decops[2], M_PI / 4, 0.0001);

  // Fully negatively correlated
  covariance.setZero();
  covariance << 4., -4., -4., 4.;
  decops = Acts::Visualization::decomposeCovariance(covariance);
  BOOST_CHECK(decops[0] == 8.);
  BOOST_CHECK(decops[1] == 0.);
  CHECK_CLOSE_ABS(decops[2], 3 * M_PI / 4, 0.0001);

  // Correlation coefficient 0.5 (off-diagonal: 3*2*0.5)
  covariance.setZero();
  covariance << 4., 2., 2., 4.;
  decops = Acts::Visualization::decomposeCovariance(covariance);
  BOOST_CHECK(decops[0] == 6.);
  BOOST_CHECK(decops[1] == 2.);
  CHECK_CLOSE_ABS(decops[2], M_PI / 4, 0.0001);

  // Correlation coefficient -0.5 & different diagonal (off-diagonal: 3*2*0.5)
  covariance.setZero();
  covariance << 9., -3., -3, 4.;
  decops = Acts::Visualization::decomposeCovariance(covariance);
  CHECK_CLOSE_ABS(decops[0], 10.4051, 0.0001);
  CHECK_CLOSE_ABS(decops[1], 2.59488, 0.0001);
  CHECK_CLOSE_ABS(decops[2], 2.70356, 0.0001);
}

/// The tests in this section are regression tests only in order
/// to catch any unexpected changes in the output format.
///
BOOST_AUTO_TEST_CASE(LineVisualizationTestsObj) {
  ObjVisualization obj;
  IVisualization::ColorType lineColor = {0, 0, 250};
  Vector3D start = {1., 1., 1.};
  Vector3D end = {4., 4., 4.};
  Acts::Visualization::drawSegment(obj, start, end, 0.1, 72, lineColor);
  obj.write("Primitives_Line");

  std::stringstream cStream;
  obj.write(cStream);
  cStream << std::setprecision(4);

  BOOST_TEST(cStream.good());
  BOOST_CHECK_EQUAL(cStream.str().size(), 4848);
}

BOOST_AUTO_TEST_CASE(ArrowVisualizationTests) {
  ObjVisualization obj;

  Vector3D start = {1., 0., 0.};
  Vector3D end = {4., 0., 0.};
  Acts::Visualization::drawArrowForward(obj, start, end, 0.1, 0.1, 3., 72,
                                        {0, 75, 0});

  start = {1., 2., 0.};
  end = {4., 2., 0.};
  Acts::Visualization::drawArrowBackward(obj, start, end, 0.1, 0.1, 3., 72,
                                         {0, 150, 0});

  start = {1., 4., 0.};
  end = {4., 4., 0.};
  Acts::Visualization::drawArrowsBoth(obj, start, end, 0.1, 0.1, 3., 72,
                                      {0, 250, 0});
  obj.write("Primitives_Arrows");

  std::stringstream cStream;
  obj.write(cStream);
  BOOST_TEST(cStream.good());
  BOOST_TEST(cStream.str().size() == 44249);
}

/// The tests in this section are regression tests only in order
/// to catch any unexpected changes in the output format.
///
BOOST_AUTO_TEST_CASE(LocalErrorVisualizationTestsObj) {
  IVisualization::ColorType surfaceColor = {120, 250, 0};

  ObjVisualization obj;
  Acts::Visualization::drawSurface(obj, *plane, gctx, Transform3D::Identity(),
                                   72, false, surfaceColor);

  // Error visualization
  IVisualization::ColorType errorColor = {250, 0, 0};

  ActsSymMatrixD<2> cov = ActsSymMatrixD<2>::Identity();
  double s0 = 0.45;
  double s1 = 1.99;
  double r01 = 0.78;
  cov << s0 * s0, r01 * s0 * s1, r01 * s0 * s1, s1 * s1;

  Vector2D lcentered{0., 0.};

  Acts::Visualization::drawCovarianceCartesian(
      obj, lcentered, cov, plane->transform(gctx), {3}, 10., 72, errorColor);

  obj.write("Primitives_CartesianError");

  std::stringstream cStream;
  cStream << std::setprecision(4);
  obj.write(cStream);
  BOOST_TEST(cStream.good());
  BOOST_CHECK_EQUAL(cStream.str().size(), 1841);
}

/// The tests in this section are regression tests only in order
/// to catch any unexpected changes in the output format.
///
BOOST_AUTO_TEST_CASE(AngularErrorVisualizationTestsObj) {
  IVisualization::ColorType surfaceColor = {120, 250, 0};

  ObjVisualization obj;

  Acts::Visualization::drawSurface(obj, *plane, gctx, Transform3D::Identity(),
                                   72, false, surfaceColor);

  // Error visualization
  IVisualization::ColorType errorColor = {250, 0, 0};

  ActsSymMatrixD<2> cov = ActsSymMatrixD<2>::Identity();
  double s0 = 0.08;
  double s1 = 0.01;
  double r01 = 0.3;
  cov << s0 * s0, r01 * s0 * s1, r01 * s0 * s1, s1 * s1;

  Vector3D origin{0., 0., 0.};
  Vector3D direction = Vector3D(1., 3., 15.).normalized();

  double directionScale = 5.;

  Acts::Visualization::drawCovarianceAngular(
      obj, origin, direction, cov, {3}, directionScale, 10., 72, errorColor);

  Acts::Visualization::drawArrowForward(
      obj, origin + 0.5 * directionScale * direction,
      origin + 1.2 * directionScale * direction, 0.02, 0.1, 5., 72,
      {10, 10, 10});

  obj.write("Primitives_AngularError");

  std::stringstream cStream;
  cStream << std::setprecision(4);
  obj.write(cStream);
  BOOST_TEST(cStream.good());
  BOOST_CHECK_EQUAL(cStream.str().size(), 15574);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts