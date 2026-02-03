// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/IndexGrid.hpp"
#include "Acts/Geometry/ReferenceGenerators.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/IndexGridNavigationPolicy.hpp"
#include "Acts/Navigation/NavigationDelegate.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <numbers>

namespace ActsTests {

using namespace Acts;

auto tContext = GeometryContext::dangerouslyDefaultConstruct();

auto tLogger = getDefaultLogger("IndexGridNavigation", Logging::VERBOSE);

BOOST_AUTO_TEST_SUITE(NavigationSuite)

BOOST_AUTO_TEST_CASE(RegularPlaneIndexGridTests) {
  // Test here:
  // - we create a grid -25, -15, -5, 5, 15, 25 in X and Y
  // - we create a surface at Z=0 spanning -4 to 4 in X and -6 to 6 in Y
  //
  // Setup the surface should thus be (without any expansion):
  // - central bin [3,3] for center reference at (0,0,0)
  // - bins [3,2], [3,3], [3,4] for polyhedron reference
  //
  // All can be modified by bin expansion or reference expansion

  // Let's create a simple plane in the XY plane
  auto planeSurface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(4., 6.));

  // x-y Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisX(-35, 35, 7);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisY(-25., 25., 5);
  Grid gridXY(Type<std::vector<std::size_t>>, std::move(axisX),
              std::move(axisY));

  TrackingVolume tVolume(Transform3::Identity(),
                         std::make_shared<CuboidVolumeBounds>(20., 20., 5.),
                         "CuboidVolume");

  tVolume.addSurface(planeSurface);

  // Indexed Surface grid
  IndexGrid<decltype(gridXY)> indexedGridXY(
      std::move(gridXY), {AxisDirection::AxisX, AxisDirection::AxisY});

  // Create a tracking volume and add the surface

  // (1a) Test with Center reference generator - no bin expansion
  IndexGridNavigationConfig centerConfig;
  centerConfig.referenceGenerator =
      std::make_shared<CenterReferenceGenerator>();

  IndexGridNavigationPolicy<decltype(gridXY)> centerNavigationPolicy(
      tContext, tVolume, *tLogger, centerConfig, indexedGridXY);

  NavigationDelegate delegate;
  BOOST_CHECK_NO_THROW(centerNavigationPolicy.connect(delegate));

  // Now initialize candidates at position (0,0,0)
  NavigationArguments navArgs;
  NavigationStream nStream;
  AppendOnlyNavigationStream navStream{nStream};

  // Address central posision
  navArgs.position = Vector3(0., 0., 0.);
  navArgs.direction = Vector3(0., 0., 1.);
  centerNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                              *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);

  // The off-central position - should yield no candidates
  nStream.reset();
  navArgs.position = Vector3(11., 11., 0.);
  navArgs.direction = Vector3(0., 0., 1.);
  centerNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                              *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 0);

  // (1b) Test with Center reference generator - with bin expansion (1,0)
  IndexGridNavigationConfig expandedConfig;
  expandedConfig.referenceGenerator =
      std::make_shared<CenterReferenceGenerator>();
  expandedConfig.binExpansion = {1u, 0u};

  IndexGridNavigationPolicy<decltype(gridXY)> expandedNavigationPolicy(
      tContext, tVolume, *tLogger, expandedConfig, indexedGridXY);

  nStream.reset();

  // The bins are expanded in X - should yield a candidate
  navArgs.position = Vector3(11., 0., 0.);
  navArgs.direction = Vector3(0., 0., 1.);
  expandedNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                                *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
  // They are not expanded in Y - should yield no candidate
  nStream.reset();
  navArgs.position = Vector3(0., 11., 0.);
  navArgs.direction = Vector3(0., 0., 1.);
  expandedNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                                *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 0);

  // (2a) Test with Polyhedron reference generator - no bin expansion
  IndexGridNavigationConfig polyConfig;
  polyConfig.referenceGenerator =
      std::make_shared<PolyhedronReferenceGenerator>();

  IndexGridNavigationPolicy<decltype(gridXY)> polyNavigationPolicy(
      tContext, tVolume, *tLogger, polyConfig, indexedGridXY);
  nStream.reset();
  // Address central posision
  navArgs.position = Vector3(0., 0., 0.);
  polyNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                            *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
  // Through the polyhedron also the bins in y before/after are filled
  nStream.reset();
  navArgs.position = Vector3(0., -7., 0.);
  polyNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                            *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
  nStream.reset();
  navArgs.position = Vector3(0., 7., 0.);
  polyNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                            *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
  // However, the bins in x before/after are not filled
  nStream.reset();
  navArgs.position = Vector3(-7., 0., 0.);
  polyNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                            *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 0);
  nStream.reset();
  navArgs.position = Vector3(7., 0., 0.);
  polyNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                            *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 0);

  // (2b) Test with Polyhedron reference generator - with reference expansion
  // (12., 0.)
  IndexGridNavigationConfig polyExpandedConfig;
  polyExpandedConfig.referenceGenerator =
      std::make_shared<PolyhedronReferenceGenerator>();
  polyExpandedConfig.referenceExpansion = {12., 0.};

  IndexGridNavigationPolicy<decltype(gridXY)> polyExpandedNavigationPolicy(
      tContext, tVolume, *tLogger, polyExpandedConfig, indexedGridXY);
  nStream.reset();
  // Address central posision - should still work
  navArgs.position = Vector3(0., 0., 0.);
  polyExpandedNavigationPolicy.initializeCandidates(tContext, navArgs,
                                                    navStream, *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
  // Sunndenly all x bins +2/-2 are filled
  nStream.reset();
  navArgs.position = Vector3(-20., 0., 0.);
  polyExpandedNavigationPolicy.initializeCandidates(tContext, navArgs,
                                                    navStream, *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
  nStream.reset();
  navArgs.position = Vector3(20., 0., 0.);
  polyExpandedNavigationPolicy.initializeCandidates(tContext, navArgs,
                                                    navStream, *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
  // While the first bin is still out of reach
  nStream.reset();
  navArgs.position = Vector3(-30., 0., 0.);
  polyExpandedNavigationPolicy.initializeCandidates(tContext, navArgs,
                                                    navStream, *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 0);
  nStream.reset();
  navArgs.position = Vector3(30., 0., 0.);
  polyExpandedNavigationPolicy.initializeCandidates(tContext, navArgs,
                                                    navStream, *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 0);
}

BOOST_AUTO_TEST_CASE(RegularCylinderIndexGridTests) {
  // Most of the tests are covered by the plane grid, we can concentrate here on
  // the phi periodicity & the projected reference generator

  // Test setup:
  // - We create a grid in z: -20 to 20 in 20 bins in z
  // - in phi: -pi to pi in 10 bins
  // - We create a surface at phi = -PI at radius R=10, spanning -2 to 2 in z

  // We create a plane surface at phi = -PI and Radius R
  double cylinderRadius = 10.;
  Vector3 surfaceCenter(-cylinderRadius, 0., 0.);

  // Local z axis is the negative x axis
  Vector3 surfaceLocalZ(-1, 0., 0.);
  // Local x axis is the global y axis
  Vector3 surfaceLocalY(0., 0., -1.);
  // Local x axis is then the global
  Vector3 surfaceLocalX(0., 1., 0);
  // Create the RotationMatrix
  RotationMatrix3 surfaceRotation;
  surfaceRotation.col(0) = surfaceLocalX;
  surfaceRotation.col(1) = surfaceLocalY;
  surfaceRotation.col(2) = surfaceLocalZ;
  // Get the surfaceTransform
  auto surfaceTransform =
      Transform3(Translation3(surfaceCenter) * surfaceRotation);

  auto planeSurface = Surface::makeShared<PlaneSurface>(
      surfaceTransform, std::make_shared<RectangleBounds>(2., 3.));

  // z-phi axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(
      -std::numbers::pi, std::numbers::pi, 10);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisZ(-20, 20, 20);
  Grid gridPhiZ(Type<std::vector<std::size_t>>, std::move(axisPhi),
                std::move(axisZ));

  TrackingVolume tVolume(Transform3::Identity(),
                         std::make_shared<CylinderVolumeBounds>(0., 15., 22.),
                         "CylinderVolume");
  tVolume.addSurface(planeSurface);

  // (1a) - Index grid with polyhedron reference generator - no expansion
  IndexGrid<decltype(gridPhiZ)> indexedgridPhiZ(
      std::move(gridPhiZ), {AxisDirection::AxisPhi, AxisDirection::AxisZ});
  IndexGridNavigationConfig polyConfig;
  polyConfig.referenceGenerator =
      std::make_shared<PolyhedronReferenceGenerator>();
  IndexGridNavigationPolicy<decltype(gridPhiZ)> polyNavigationPolicy(
      tContext, tVolume, *tLogger, polyConfig, indexedgridPhiZ);

  NavigationDelegate delegate;
  BOOST_CHECK_NO_THROW(polyNavigationPolicy.connect(delegate));

  // Now initialize candidates at position (R,0,0) - opposite should. not lead
  // to a candidate
  NavigationArguments navArgs;
  NavigationStream nStream;
  AppendOnlyNavigationStream navStream{nStream};
  navArgs.position = Vector3(cylinderRadius, 0., 0.);
  navArgs.direction = Vector3(1., 0., 0.);
  polyNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                            *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 0);
  // However, a position at (-R,0,0) should yield a candidate
  nStream.reset();
  navArgs.position = Vector3(-cylinderRadius, 0., 0.);
  navArgs.direction = Vector3(-1., 0., 0.);
  polyNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                            *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
  // A candidate a off in phi - no results
  nStream.reset();
  navArgs.position =
      Vector3(cylinderRadius * std::cos(-std::numbers::pi + 0.8),
              cylinderRadius * std::sin(-std::numbers::pi + 0.8), 0.);
  polyNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                            *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 0);

  // A candidate a off in z range - no result
  nStream.reset();
  navArgs.position = Vector3(-cylinderRadius, 0., 4.);
  polyNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                            *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 0);

  // Let's test a reference surface projection generator
  // Create a reference surface with bigger radius, should expand z but not phi
  double referenceSurfaceRadius = 15.;
  auto referenceCylinder = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(),
      std::make_shared<CylinderBounds>(referenceSurfaceRadius, 22.));

  IndexGridNavigationConfig projectedConfig;
  auto projectedReferenceGenerator =
      std::make_shared<ProjectedReferenceGenerator>();
  projectedReferenceGenerator->nSegements = 1;
  projectedReferenceGenerator->expansionValue = 0.0;
  projectedReferenceGenerator->referenceSurface = referenceCylinder;
  projectedReferenceGenerator->luminousRegion = {Vector3(0., 0., 0.)};

  projectedConfig.referenceGenerator = projectedReferenceGenerator;
  IndexGridNavigationPolicy<decltype(gridPhiZ)> projectedNavigationPolicy(
      tContext, tVolume, *tLogger, projectedConfig, indexedgridPhiZ);

  // A candidate a off in phi - still outside
  nStream.reset();
  navArgs.position =
      Vector3(cylinderRadius * std::cos(-std::numbers::pi + 0.8),
              cylinderRadius * std::sin(-std::numbers::pi + 0.8), 0.);
  projectedNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                                 *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 0);

  // A candidate a off in z range - should now yield a candidate now
  nStream.reset();
  navArgs.position = Vector3(-cylinderRadius, 0., 4.);
  projectedNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                                 *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);

  // Now reference generator that also carries the reference surface
  IndexGridNavigationConfig projectedWithSurfaceConfig;
  projectedWithSurfaceConfig.surface = referenceCylinder;
  auto projectedWithSurfaceReferenceGenerator =
      std::make_shared<ProjectedReferenceGenerator>();
  projectedWithSurfaceReferenceGenerator->nSegements = 1;
  projectedWithSurfaceReferenceGenerator->expansionValue = 0.0;
  projectedWithSurfaceReferenceGenerator->referenceSurface = referenceCylinder;
  projectedWithSurfaceReferenceGenerator->luminousRegion = {
      Vector3(0., 0., 0.)};
  projectedWithSurfaceConfig.referenceGenerator =
      projectedWithSurfaceReferenceGenerator;
  IndexGridNavigationPolicy<decltype(gridPhiZ)>
      projectedWithSurfaceNavigationPolicy(tContext, tVolume, *tLogger,
                                           projectedWithSurfaceConfig,
                                           indexedgridPhiZ);
  // A candidate a off in phi - we get the phi surface now as well
  nStream.reset();
  navArgs.position =
      Vector3(cylinderRadius * std::cos(-std::numbers::pi + 0.8),
              cylinderRadius * std::sin(-std::numbers::pi + 0.8), 0.);
  projectedWithSurfaceNavigationPolicy.initializeCandidates(
      tContext, navArgs, navStream, *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
}

BOOST_AUTO_TEST_CASE(RegularDiscIndexGridTests) {
  // Most of the test cases are covered by the plane and cylinder grid and
  // planar grid we can concentrate here on multiple fillings

  // Test setup:
  // - We create a grid in r: 5 to 25 in 4 bins
  // - in phi: -pi to pi in 10 bins
  // - We create a few disc sectoral surfaces and let them overlap in r and phi

  // This surface spans in some bins around phi = 0
  auto surface0 = Surface::makeShared<DiscSurface>(
      Transform3::Identity(),
      std::make_shared<RadialBounds>(5., 12., std::numbers::pi / 5, 0.));

  // Tiny surface only in the neighboring bin
  auto surface1 = Surface::makeShared<DiscSurface>(
      Transform3::Identity(),
      std::make_shared<RadialBounds>(6., 8., std::numbers::pi / 12,
                                     1.5 * std::numbers::pi / 5));

  double innerRadius = 5.;
  double outerRadius = 25.;

  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisR(innerRadius,
                                                             outerRadius, 4);
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(
      -std::numbers::pi, std::numbers::pi, 10);
  Grid gridRPhi(Type<std::vector<std::size_t>>, std::move(axisR),
                std::move(axisPhi));

  auto tVolume = TrackingVolume(Transform3::Identity(),
                                std::make_shared<CylinderVolumeBounds>(
                                    innerRadius - 1., outerRadius + 1., 5.),
                                "DiscVolume");
  tVolume.addSurface(surface0);
  tVolume.addSurface(surface1);

  // (1a) - Index grid with center reference generator - no expansion
  IndexGrid<decltype(gridRPhi)> indexedGridRPhi(
      std::move(gridRPhi), {AxisDirection::AxisR, AxisDirection::AxisPhi});
  IndexGridNavigationConfig centerConfig;
  centerConfig.referenceGenerator =
      std::make_shared<PolyhedronReferenceGenerator>();
  IndexGridNavigationPolicy<decltype(gridRPhi)> centerNavigationPolicy(
      tContext, tVolume, *tLogger, centerConfig, indexedGridRPhi);

  NavigationDelegate delegate;
  BOOST_CHECK_NO_THROW(centerNavigationPolicy.connect(delegate));

  // Now initialize candidates at position (R,0,0) - should yield only surface0
  NavigationArguments navArgs;
  NavigationStream nStream;
  AppendOnlyNavigationStream navStream{nStream};
  navArgs.position = Vector3(7.5, 0., 0.);
  navArgs.direction = Vector3(1., 0., 0.);
  centerNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                              *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
  BOOST_CHECK(&nStream.currentCandidate().surface() == surface0.get());
  // Check the neighboring bin, which should not be filled by surface0 and
  // surface1
  nStream.reset();
  navArgs.position = Vector3(7.5 * std::cos(1.1 * std::numbers::pi / 5),
                             7.5 * std::sin(1.1 * std::numbers::pi / 5), 0.);
  centerNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                              *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 2);
}

BOOST_AUTO_TEST_CASE(RegularRingIndexGridTests) {
  // This tests the one dimensional closed grid in phi only
  // Test setup:
  // - We create a grid in phi: -pi to pi in 10 bins

  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(
      -std::numbers::pi, std::numbers::pi, 10);
  Grid gridPhi(Type<std::vector<std::size_t>>, std::move(axisPhi));

  auto surface = Surface::makeShared<DiscSurface>(
      Transform3::Identity(),
      std::make_shared<RadialBounds>(6., 14., std::numbers::pi / 5, 0.));

  auto tVolume = TrackingVolume(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(0., 15., 5.), "RingVolume");
  tVolume.addSurface(surface);
  // (1a) - Index grid with center reference generator - no expansion
  IndexGrid<decltype(gridPhi)> indexedGridPhi(std::move(gridPhi),
                                              {AxisDirection::AxisPhi});
  IndexGridNavigationConfig centerConfig;
  centerConfig.binExpansion = {0u, 0u};
  centerConfig.referenceGenerator =
      std::make_shared<CenterReferenceGenerator>();
  BOOST_CHECK_THROW(
      {
        IndexGridNavigationPolicy<decltype(gridPhi)>
            centerNavigationPolicyThrow(tContext, tVolume, *tLogger,
                                        centerConfig, indexedGridPhi);
      },
      std::runtime_error);

  centerConfig.binExpansion = {1u};

  IndexGridNavigationPolicy<decltype(gridPhi)> centerNavigationPolicy(
      tContext, tVolume, *tLogger, centerConfig, indexedGridPhi);

  NavigationDelegate delegate;
  BOOST_CHECK_NO_THROW(centerNavigationPolicy.connect(delegate));

  // Now initialize candidates at position (R,0,0) - should yield the surface
  NavigationArguments navArgs;
  NavigationStream nStream;
  AppendOnlyNavigationStream navStream{nStream};
  navArgs.position = Vector3(10., 0., 0.);
  navArgs.direction = Vector3(1., 0., 0.);
  centerNavigationPolicy.initializeCandidates(tContext, navArgs, navStream,
                                              *tLogger);
  BOOST_CHECK_EQUAL(nStream.candidates().size(), 1);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
