// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/SurfaceAxisResolution.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/AxisFactory.hpp"
#include "Acts/Utilities/MultiAxisFactory.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <numbers>
#include <stdexcept>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(SurfacesSuite)

BOOST_AUTO_TEST_CASE(SurfaceAxisResolutionCylinder) {
  using enum AxisDirection;
  using enum AxisBoundaryType;

  // Full cylinder: phi-like axes are closed
  auto fullCylinder = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(25., 50.));
  std::array<AxisDirection, 2> expected = {AxisRPhi, AxisZ};
  const Surface& cylinderRef = *fullCylinder;
  BOOST_CHECK(cylinderRef.localAxes() == expected);

  auto rphi = surfaceAxisResolution(*fullCylinder, AxisRPhi);
  CHECK_CLOSE_ABS(rphi.min, -25. * std::numbers::pi, 1e-12);
  CHECK_CLOSE_ABS(rphi.max, 25. * std::numbers::pi, 1e-12);
  BOOST_CHECK_EQUAL(rphi.boundaryType, Closed);

  auto phi = surfaceAxisResolution(*fullCylinder, AxisPhi);
  CHECK_CLOSE_ABS(phi.min, -std::numbers::pi, 1e-12);
  BOOST_CHECK_EQUAL(phi.boundaryType, Closed);

  auto z = surfaceAxisResolution(*fullCylinder, AxisZ);
  CHECK_CLOSE_ABS(z.min, -50., 1e-12);
  CHECK_CLOSE_ABS(z.max, 50., 1e-12);
  BOOST_CHECK_EQUAL(z.boundaryType, Bound);

  BOOST_CHECK_THROW(surfaceAxisResolution(*fullCylinder, AxisR),
                    std::invalid_argument);

  // Sectoral cylinder: phi-like axes are bound
  auto sector = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(),
      std::make_shared<CylinderBounds>(25., 50., std::numbers::pi / 4.));
  auto sectorRPhi = surfaceAxisResolution(*sector, AxisRPhi);
  CHECK_CLOSE_ABS(sectorRPhi.min, -25. * std::numbers::pi / 4., 1e-12);
  BOOST_CHECK_EQUAL(sectorRPhi.boundaryType, Bound);
}

BOOST_AUTO_TEST_CASE(SurfaceAxisResolutionDiscAndPlane) {
  using enum AxisDirection;
  using enum AxisBoundaryType;

  auto disc = Surface::makeShared<DiscSurface>(
      Transform3::Identity(), std::make_shared<RadialBounds>(50., 75.));
  std::array<AxisDirection, 2> expectedDisc = {AxisR, AxisPhi};
  const Surface& discRef = *disc;
  BOOST_CHECK(discRef.localAxes() == expectedDisc);
  auto r = surfaceAxisResolution(*disc, AxisR);
  CHECK_CLOSE_ABS(r.min, 50., 1e-12);
  CHECK_CLOSE_ABS(r.max, 75., 1e-12);
  BOOST_CHECK_EQUAL(r.boundaryType, Bound);
  BOOST_CHECK_EQUAL(surfaceAxisResolution(*disc, AxisPhi).boundaryType, Closed);

  auto rectangle = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(20., 30.));
  std::array<AxisDirection, 2> expectedPlane = {AxisX, AxisY};
  const Surface& rectangleRef = *rectangle;
  BOOST_CHECK(rectangleRef.localAxes() == expectedPlane);
  auto x = surfaceAxisResolution(*rectangle, AxisX);
  CHECK_CLOSE_ABS(x.min, -20., 1e-12);
  CHECK_CLOSE_ABS(x.max, 20., 1e-12);
  BOOST_CHECK_EQUAL(x.boundaryType, Bound);

  auto trapezoid = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<TrapezoidBounds>(5., 15., 30.));
  auto tx = surfaceAxisResolution(*trapezoid, AxisX);
  CHECK_CLOSE_ABS(tx.min, -15., 1e-12);
  CHECK_CLOSE_ABS(tx.max, 15., 1e-12);
  auto ty = surfaceAxisResolution(*trapezoid, AxisY);
  CHECK_CLOSE_ABS(ty.max, 30., 1e-12);

  // Unsupported surface
  auto straw =
      Surface::makeShared<StrawSurface>(Transform3::Identity(), 5., 100.);
  BOOST_CHECK_THROW(surfaceAxisResolution(*straw, AxisZ),
                    std::invalid_argument);
  MultiAxisFactory strawBinning({AxisFactory::DeferredEquidistant(4)});
  BOOST_CHECK_THROW(resolveAxes(strawBinning, *straw), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ResolveAxesAgainstSurface) {
  using enum AxisDirection;
  using enum AxisBoundaryType;

  auto cylinder = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(30., 100.));

  // Positional matching without directions
  MultiAxisFactory positional({AxisFactory::DeferredEquidistant(10),
                               AxisFactory::DeferredEquidistant(20)});
  auto axes = resolveAxes(positional, *cylinder);
  BOOST_CHECK_EQUAL(axes.size(), 2);
  BOOST_CHECK(axes[0]->getDirection() == AxisRPhi);
  BOOST_CHECK_EQUAL(axes[0]->getBoundaryType(), Closed);
  BOOST_CHECK_EQUAL(axes[0]->getNBins(), 10);
  BOOST_CHECK(axes[1]->getDirection() == AxisZ);
  BOOST_CHECK_EQUAL(axes[1]->getBoundaryType(), Bound);
  CHECK_CLOSE_ABS(axes[1]->getMin(), -100., 1e-12);

  // Directed matching in canonical order
  MultiAxisFactory directed({AxisFactory::DeferredEquidistant(10, AxisRPhi),
                             AxisFactory::DeferredEquidistant(20, AxisZ)});
  auto directedAxes = resolveAxes(directed, *cylinder);
  BOOST_CHECK_EQUAL(directedAxes[0]->getNBins(), 10);
  BOOST_CHECK_EQUAL(directedAxes[1]->getNBins(), 20);

  // Swapped order is re-ordered to match the surface
  MultiAxisFactory swapped({AxisFactory::DeferredEquidistant(20, AxisZ),
                            AxisFactory::DeferredEquidistant(10, AxisRPhi)});
  auto reordered = resolveAxes(swapped, *cylinder);
  BOOST_CHECK(reordered[0]->getDirection() == AxisRPhi);
  BOOST_CHECK_EQUAL(reordered[0]->getNBins(), 10);
  BOOST_CHECK(reordered[1]->getDirection() == AxisZ);
  BOOST_CHECK_EQUAL(reordered[1]->getNBins(), 20);

  // 1D binning picks its direction
  MultiAxisFactory oneD({AxisFactory::DeferredEquidistant(8, AxisZ)});
  auto oneDAxes = resolveAxes(oneD, *cylinder);
  BOOST_CHECK_EQUAL(oneDAxes.size(), 1);
  BOOST_CHECK(oneDAxes[0]->getDirection() == AxisZ);

  // 1D binning without a direction binds to the first canonical direction
  MultiAxisFactory oneDFree({AxisFactory::DeferredEquidistant(8)});
  auto oneDFreeAxes = resolveAxes(oneDFree, *cylinder);
  BOOST_CHECK(oneDFreeAxes[0]->getDirection() == AxisRPhi);

  // Deferred variable binning is scaled onto the surface range
  MultiAxisFactory variable(
      {AxisFactory::DeferredVariable({0., 0.25, 1.}, AxisZ)});
  auto variableAxes = resolveAxes(variable, *cylinder);
  auto edges = variableAxes[0]->getBinEdges();
  CHECK_CLOSE_ABS(edges[0], -100., 1e-12);
  CHECK_CLOSE_ABS(edges[1], -50., 1e-12);
  CHECK_CLOSE_ABS(edges[2], 100., 1e-12);

  // Mismatching directions are rejected
  MultiAxisFactory wrongDirs({AxisFactory::DeferredEquidistant(10, AxisR),
                              AxisFactory::DeferredEquidistant(20, AxisPhi)});
  BOOST_CHECK_THROW(resolveAxes(wrongDirs, *cylinder), std::invalid_argument);

  // Duplicate directions are rejected
  MultiAxisFactory duplicate({AxisFactory::DeferredEquidistant(10, AxisZ),
                              AxisFactory::DeferredEquidistant(20, AxisZ)});
  BOOST_CHECK_THROW(resolveAxes(duplicate, *cylinder), std::invalid_argument);

  // Mixed directed and undirected axes are rejected
  MultiAxisFactory mixed({AxisFactory::DeferredEquidistant(10, AxisRPhi),
                          AxisFactory::DeferredEquidistant(20)});
  BOOST_CHECK_THROW(resolveAxes(mixed, *cylinder), std::invalid_argument);

  // Too many dimensions are rejected
  MultiAxisFactory tooMany({AxisFactory::DeferredEquidistant(2),
                            AxisFactory::DeferredEquidistant(2),
                            AxisFactory::DeferredEquidistant(2)});
  BOOST_CHECK_THROW(resolveAxes(tooMany, *cylinder), std::invalid_argument);

  // Fully specified descriptions are rejected (fail fast)
  MultiAxisFactory full({AxisFactory::Equidistant(Bound, 0., 1., 10),
                         AxisFactory::DeferredEquidistant(20)});
  BOOST_CHECK_THROW(resolveAxes(full, *cylinder), std::domain_error);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
