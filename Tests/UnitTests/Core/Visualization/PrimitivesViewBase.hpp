// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Visualization/EventDataView.hpp"
#include "Acts/Visualization/GeometryView.hpp"
#include "Acts/Visualization/IVisualization.hpp"
#include "Acts/Visualization/ObjVisualization.hpp"
#include "Acts/Visualization/PlyVisualization.hpp"

#include <fstream>
#include <sstream>
#include <string>

namespace Acts {

namespace PrimitivesViewTest {

// Test on a plane
auto identity = std::make_shared<Transform3D>(Transform3D::Identity());
auto rectangle = std::make_shared<RectangleBounds>(10., 10.);
auto plane = Surface::makeShared<PlaneSurface>(identity, rectangle);

// Test context
GeometryContext gctx = GeometryContext();

/// Helper method to visualiza all types of surfaces
///
/// @param helper The visualziation helper
///
/// @return an overall string including all written output
static inline std::string run(IVisualization& helper) {
  std::stringstream ss;

  ViewConfig lineView({0, 0, 255});
  lineView.lineThickness = 0.1;

  // Line visualization ------------------------------------------------
  Vector3D start = {1., 1., 1.};
  Vector3D end = {4., 4., 4.};
  Acts::GeometryView::drawSegment(helper, start, end);
  helper.write("Primitives_Line");
  helper.write(ss);
  helper.clear();

  // Arrows visualization ------------------------------------------------
  start = {1., 0., 0.};
  end = {4., 0., 0.};
  Acts::GeometryView::drawArrowForward(helper, start, end, 3., 2., lineView);

  start = {1., 2., 0.};
  end = {4., 2., 0.};
  Acts::GeometryView::drawArrowBackward(helper, start, end, 3., 2., lineView);

  start = {1., 4., 0.};
  end = {4., 4., 0.};
  Acts::GeometryView::drawArrowsBoth(helper, start, end, 3., 2., lineView);

  helper.write("Primitives_Arrows");
  helper.write(ss);
  helper.clear();

  // Error visualization: local ---------------------------------------------
  Acts::GeometryView::drawSurface(helper, *plane, gctx);

  ViewConfig errorVis({250, 0, 0});
  errorVis.lineThickness = 0.025;

  ActsSymMatrixD<2> cov = ActsSymMatrixD<2>::Identity();
  double s0 = 0.75;
  double s1 = 1.99;
  double r01 = 0.78;
  cov << s0 * s0, r01 * s0 * s1, r01 * s0 * s1, s1 * s1;

  Vector2D lcentered{0., 0.};
  Acts::EventDataView::drawCovarianceCartesian(
      helper, lcentered, cov, plane->transform(gctx), 1.0, errorVis);

  helper.write("Primitives_CartesianError");
  helper.write(ss);
  helper.clear();

  // Error visualization: angular ---------------------------------------------
  Acts::GeometryView::drawSurface(helper, *plane, gctx);
  cov = ActsSymMatrixD<2>::Identity();
  s0 = 0.08;
  s1 = 0.035;
  r01 = 0.3;
  cov << s0 * s0, r01 * s0 * s1, r01 * s0 * s1, s1 * s1;

  Vector3D origin{0., 0., 0.};
  Vector3D direction = Vector3D(1., 3., 15.).normalized();

  double directionScale = 5.;

  Acts::EventDataView::drawCovarianceAngular(helper, origin, direction, cov,
                                             directionScale, 10., errorVis);

  Acts::GeometryView::drawArrowForward(
      helper, origin + 0.5 * directionScale * direction,
      origin + 1.2 * directionScale * direction, 3., 2., errorVis);

  helper.write("Primitives_AngularError");
  helper.write(ss);
  helper.clear();

  return ss.str();
}

}  // namespace PrimitivesViewTest
}  // namespace Acts