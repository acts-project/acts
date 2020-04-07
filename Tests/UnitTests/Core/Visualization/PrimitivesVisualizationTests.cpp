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
#include "Acts/Visualization/EventDataVisualization.hpp"
#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/IVisualization.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"
#include "Acts/Visualization/PlyVisualization.hpp"

#include <fstream>
#include <iostream>

namespace Acts {

namespace Test {

BOOST_AUTO_TEST_SUITE(VisualizationTests)

// Test on a plane
auto identity = std::make_shared<Transform3D>(Transform3D::Identity());
auto rectangle = std::make_shared<RectangleBounds>(10., 10.);
auto plane = Surface::makeShared<PlaneSurface>(identity, rectangle);

// Test context
GeometryContext gctx = GeometryContext();

BOOST_AUTO_TEST_CASE(LineVisualizationTestsObj) {
  ObjVisualization obj;
  IVisualization::ColorType lineColor = {0, 0, 250};
  Vector3D start = {1., 1., 1.};
  Vector3D end = {4., 4., 4.};
  Visualization::drawSegment(obj, start, end, 0.1, 72, lineColor);
  obj.write("Primitives_Line");
}

BOOST_AUTO_TEST_CASE(ArrowVisualizationTests) {
  ObjVisualization obj;

  Vector3D start = {1., 0., 0.};
  Vector3D end = {4., 0., 0.};
  Visualization::drawArrowForward(obj, start, end, 0.1, 0.1, 3., 72,
                                  {0, 75, 0});

  start = {1., 2., 0.};
  end = {4., 2., 0.};
  Visualization::drawArrowBackward(obj, start, end, 0.1, 0.1, 3., 72,
                                   {0, 150, 0});

  start = {1., 4., 0.};
  end = {4., 4., 0.};
  Visualization::drawArrowsBoth(obj, start, end, 0.1, 0.1, 3., 72, {0, 250, 0});
  obj.write("Primitives_Arrows");
}

BOOST_AUTO_TEST_CASE(LocalErrorVisualizationTestsObj) {
  IVisualization::ColorType surfaceColor = {120, 250, 0};

  ObjVisualization obj;

  Visualization::drawSurface(obj, *plane, gctx, Transform3D::Identity(), 72,
                             false, surfaceColor);

  // Error visualization
  IVisualization::ColorType errorColor = {250, 0, 0};

  ActsSymMatrixD<2> cov = ActsSymMatrixD<2>::Identity();
  double s0 = 0.45;
  double s1 = 1.99;
  double r01 = 0.78;
  cov << s0 * s0, r01 * s0 * s1, r01 * s0 * s1, s1 * s1;

  Vector2D lcentered{0., 0.};

  Visualization::drawCovarianceCartesian(
      obj, lcentered, cov, plane->transform(gctx), {3}, 10., 72, errorColor);

  obj.write("Primitives_CartesianError");
}

BOOST_AUTO_TEST_CASE(AngularErrorVisualizationTestsObj) {
  IVisualization::ColorType surfaceColor = {120, 250, 0};

  ObjVisualization obj;

  Visualization::drawSurface(obj, *plane, gctx, Transform3D::Identity(), 72,
                             false, surfaceColor);

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

  Visualization::drawCovarianceAngular(obj, origin, direction, cov, {3},
                                       directionScale, 10., 72, errorColor);

  Visualization::drawArrowForward(obj,
                                  origin + 0.5 * directionScale * direction,
                                  origin + 1.2 * directionScale * direction,
                                  0.02, 0.1, 5., 72, {10, 10, 10});

  obj.write("Primitives_AngularError");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts