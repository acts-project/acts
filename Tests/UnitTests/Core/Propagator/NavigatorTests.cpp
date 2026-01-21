// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <cstddef>
#include <memory>
#include <string>

namespace bdata = boost::unit_test::data;

using namespace Acts;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::perp;

namespace ActsTests {

// Create a test context
GeometryContext tgContext = GeometryContext::dangerouslyDefaultConstruct();
MagneticFieldContext mfContext = MagneticFieldContext();

void step(Vector3& pos, const Vector3& dir, double stepSize) {
  pos += stepSize * dir;
}

void step(Vector3& pos, const Vector3& dir, const Surface& surface) {
  Intersection3D intersection =
      surface.intersect(tgContext, pos, dir).closestForward();
  step(pos, dir, intersection.pathLength());
}

void step(Vector3& pos, const Vector3& dir, const NavigationTarget& target) {
  step(pos, dir, target.surface());
}

/// @brief Method for testing vectors in @c Navigator::State
///
/// @param [in] state Navigator state
/// @param [in] navSurf Number of navigation surfaces
/// @param [in] navLay Number of navigation layers
/// @param [in] navBound Number of navigation boundaries
/// @param [in] extSurf Number of external surfaces
bool testNavigatorStateVectors(Navigator::State& state, std::size_t navSurf,
                               std::size_t navLay, std::size_t navBound) {
  return ((state.navSurfaces.size() == navSurf) &&
          (state.navLayers.size() == navLay) &&
          (state.navBoundaries.size() == navBound));
}

/// @brief Method for testing pointers in @c Navigator::State
///
/// @param [in] state Navigation state
/// @param [in] startVol Start volume
/// @param [in] startLay Start layer
/// @param [in] startSurf Start surface
/// @param [in] currSurf Current surface
/// @param [in] currVol Current volume
/// @param [in] targetSurf Target surface
bool testNavigatorStatePointers(Navigator::State& state,
                                const TrackingVolume* startVol,
                                const Layer* startLay, const Surface* startSurf,
                                const Surface* currSurf,
                                const TrackingVolume* currVol,
                                const Surface* targetSurf) {
  std::cout << "startVol: " << startVol << " startLay: " << startLay
            << " startSurf: " << startSurf << " currSurf: " << currSurf
            << " currVol: " << currVol << " targetSurf: " << targetSurf
            << std::endl;

  std::cout << "state.startVolume: " << state.startVolume
            << " state.startLayer: " << state.startLayer
            << " state.startSurface: " << state.startSurface
            << " state.currentSurface: " << state.currentSurface
            << " state.currentVolume: " << state.currentVolume
            << " state.targetSurface: " << state.targetSurface << std::endl;

  return (
      (state.startVolume == startVol) && (state.startLayer == startLay) &&
      (state.startSurface == startSurf) && (state.currentSurface == currSurf) &&
      (state.currentVolume == currVol) && (state.targetSurface == targetSurf));
}

// the surface cache & the creation of the geometry
CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

const double Bz = 2_T;
auto bField = std::make_shared<ConstantBField>(Vector3{0, 0, Bz});

Logging::Level logLevel = Logging::INFO;

BOOST_AUTO_TEST_SUITE(PropagatorSuite)

BOOST_AUTO_TEST_CASE(Navigator_status_methods) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("NavigatorTest", logLevel));

  // position and direction vector
  Vector3 position = Vector3::Zero();
  Vector3 direction = Vector3(1., 1., 0).normalized();

  ACTS_INFO("(1) Test for inactivity");
  ACTS_INFO("    a) Run without anything present");
  {
    auto bounds = std::make_shared<CylinderVolumeBounds>(10, 20, 20);
    auto tvol = std::make_shared<TrackingVolume>(Transform3::Identity(), bounds,
                                                 "Undefined");

    auto tgeo = std::make_shared<TrackingGeometry>(tvol);

    Navigator::Config navCfg;
    navCfg.resolveSensitive = false;
    navCfg.resolveMaterial = false;
    navCfg.resolvePassive = false;
    navCfg.trackingGeometry = tgeo;
    Navigator navigator{navCfg};

    Navigator::Options options(tgContext);

    Navigator::State state = navigator.makeState(options);

    BOOST_CHECK(testNavigatorStateVectors(state, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state, nullptr, nullptr, nullptr,
                                           nullptr, nullptr, nullptr));
  }

  ACTS_INFO("    b) Run with geometry but without resolving");
  {
    Navigator::Config navCfg;
    navCfg.resolveSensitive = false;
    navCfg.resolveMaterial = false;
    navCfg.resolvePassive = false;
    navCfg.trackingGeometry = tGeometry;
    Navigator navigator{navCfg};

    Navigator::Options options(tgContext);

    Navigator::State state = navigator.makeState(options);

    BOOST_CHECK(testNavigatorStateVectors(state, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state, nullptr, nullptr, nullptr,
                                           nullptr, nullptr, nullptr));
  }

  ACTS_INFO(
      "    c) Run with geometry and resolving but broken navigation for "
      "various reasons");
  {
    Navigator::Config navCfg;
    navCfg.resolveSensitive = true;
    navCfg.resolveMaterial = true;
    navCfg.resolvePassive = true;
    navCfg.trackingGeometry = tGeometry;
    Navigator navigator{navCfg};

    Navigator::Options options(tgContext);

    Navigator::State state = navigator.makeState(options);

    ACTS_INFO("        i) Because target is reached");
    state.navigationBreak = true;
    navigator.nextTarget(state, position, direction);
    BOOST_CHECK(testNavigatorStateVectors(state, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state, nullptr, nullptr, nullptr,
                                           nullptr, nullptr, nullptr));

    ACTS_INFO("        ii) Because of no target surface");
    state.targetSurface = nullptr;
    navigator.nextTarget(state, position, direction);
    BOOST_CHECK(testNavigatorStateVectors(state, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state, nullptr, nullptr, nullptr,
                                           nullptr, nullptr, nullptr));

    ACTS_INFO("        iii) Because the target surface is reached");
    auto beamline = Surface::makeShared<PerigeeSurface>(Vector3::Zero());
    const Surface* startSurf = beamline.get();
    position = startSurf->center(tgContext);
    const TrackingVolume* startVol =
        tGeometry->lowestTrackingVolume(tgContext, position);
    const Layer* startLay = startVol->associatedLayer(tgContext, position);
    state.options.startSurface = startSurf;
    state.options.targetSurface = startSurf;
    BOOST_CHECK(
        navigator.initialize(state, position, direction, Direction::Forward())
            .ok());
    BOOST_CHECK(testNavigatorStateVectors(state, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state, startVol, startLay, startSurf,
                                           startSurf, startVol, startSurf));

    ACTS_INFO("(2) Test the initialisation");
    ACTS_INFO("    a) Initialise without additional information");
    state = navigator.makeState(options);
    position = Vector3::Zero();
    startVol = tGeometry->lowestTrackingVolume(tgContext, position);
    startLay = startVol->associatedLayer(tgContext, position);
    BOOST_CHECK(
        navigator.initialize(state, position, direction, Direction::Forward())
            .ok());
    BOOST_CHECK(testNavigatorStateVectors(state, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state, startVol, startLay, nullptr,
                                           nullptr, startVol, nullptr));

    ACTS_INFO("    b) Initialise having a start surface");
    state = navigator.makeState(options);
    state.options.startSurface = startSurf;
    BOOST_CHECK(
        navigator.initialize(state, position, direction, Direction::Forward())
            .ok());
    BOOST_CHECK(testNavigatorStateVectors(state, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state, startVol, startLay, startSurf,
                                           startSurf, startVol, nullptr));

    ACTS_INFO("    c) Initialise having a start volume");
    state = navigator.makeState(options);
    state.startVolume = startVol;
    BOOST_CHECK(
        navigator.initialize(state, position, direction, Direction::Forward())
            .ok());
    BOOST_CHECK(testNavigatorStateVectors(state, 0u, 0u, 0u));
    BOOST_CHECK(testNavigatorStatePointers(state, startVol, startLay, nullptr,
                                           nullptr, startVol, nullptr));
  }
}

BOOST_AUTO_TEST_CASE(Navigator_target_methods) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("NavigatorTest", logLevel));

  // create a navigator
  Navigator::Config navCfg;
  navCfg.trackingGeometry = tGeometry;
  navCfg.resolveSensitive = true;
  navCfg.resolveMaterial = true;
  navCfg.resolvePassive = false;
  Navigator navigator{navCfg};

  Navigator::Options options(tgContext);

  Navigator::State state = navigator.makeState(options);

  // position and direction vector
  Vector3 position = Vector3::Zero();
  Vector3 direction = Vector3(1., 1., 0).normalized();

  // forward navigation ----------------------------------------------
  ACTS_INFO("<<<<<<<<<<<<<<<<<<<<< FORWARD NAVIGATION >>>>>>>>>>>>>>>>>>");

  // (1) Initialization navigation from start point
  // - this will call resolveLayers() as well
  // - and thus should call a return to the stepper
  BOOST_CHECK(
      navigator.initialize(state, position, direction, Direction::Forward())
          .ok());
  // Check that the currentVolume is set
  BOOST_CHECK_NE(state.currentVolume, nullptr);
  // Check that the currentVolume is the startVolume
  BOOST_CHECK_EQUAL(state.currentVolume, state.startVolume);
  // Check that the currentSurface is reset to:
  BOOST_CHECK_EQUAL(state.currentSurface, nullptr);
  // No layer has been found
  BOOST_CHECK_EQUAL(state.navLayers.size(), 0u);

  // Estimate the next target
  NavigationTarget target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());
  // A layer has been found
  BOOST_CHECK_EQUAL(state.navLayers.size(), 1u);
  // The index should points to the begin
  BOOST_CHECK_EQUAL(state.navLayerIndex.value(), 0);
  // Check the target is correct
  BOOST_CHECK_EQUAL(&target.surface(), &state.navLayer().surface());
  // Intersect the target
  Intersection3D targetIntersection =
      target.surface()
          .intersect(tgContext, position, direction)
          .closestForward();
  // Cache the beam pipe radius
  double beamPipeR = perp(state.navLayer().position());
  // step size has been updated
  CHECK_CLOSE_ABS(targetIntersection.pathLength(), beamPipeR,
                  s_onSurfaceTolerance);

  ACTS_INFO("<<< Test 1a >>> initialize at " << toString(position));

  // Do the step towards the beam pipe
  step(position, direction, target);

  // (2) re-entering navigator:
  // POST STEP
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  // Check that the currentVolume is the still startVolume
  BOOST_CHECK_EQUAL(state.currentVolume, state.startVolume);
  // The layer number has not changed
  BOOST_CHECK_EQUAL(state.navLayers.size(), 1u);
  // The index still points to the begin
  BOOST_CHECK_EQUAL(state.navLayerIndex.value(), 0);

  ACTS_INFO("<<< Test 1b >>> step to the BeamPipe at  " << toString(position));

  // Estimate the next target
  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());

  // Do the step towards the boundary
  step(position, direction, target);

  // (3) re-entering navigator:
  // POST STEP
  navigator.handleSurfaceReached(state, position, direction, target.surface());

  ACTS_INFO("<<< Test 1c >>> step to the Boundary at  " << toString(position));

  // Estimate the next target
  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());
  // Intersect the target
  targetIntersection = target.surface()
                           .intersect(tgContext, position, direction)
                           .closestForward();

  // positive return: do the step
  step(position, direction, target);

  // (4) re-entering navigator:
  // POST STEP
  navigator.handleSurfaceReached(state, position, direction, target.surface());

  ACTS_INFO("<<< Test 1d >>> step to 1st layer at  " << toString(position));

  // Estimate the next target
  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());

  // Step through the surfaces on first layer
  for (std::size_t isf = 0; isf < 5; ++isf) {
    step(position, direction, target);
    // (5-9) re-entering navigator:
    // POST STEP
    navigator.handleSurfaceReached(state, position, direction,
                                   target.surface());
    // ACTORS - ABORTERS - PRE STEP
    target = navigator.nextTarget(state, position, direction);
    BOOST_CHECK(!target.isNone());

    ACTS_INFO("<<< Test 1e-1i >>> step within 1st layer at  "
              << toString(position));
  }

  // positive return: do the step
  step(position, direction, target);
  // (10) re-entering navigator:
  // POST STEP
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  // ACTORS - ABORTERS - PRE STEP
  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());

  ACTS_INFO("<<< Test 1j >>> step to 2nd layer at  " << toString(position));

  // Step through the surfaces on second layer
  for (std::size_t isf = 0; isf < 5; ++isf) {
    step(position, direction, target);
    // (11-15) re-entering navigator:
    // POST STEP
    navigator.handleSurfaceReached(state, position, direction,
                                   target.surface());
    // ACTORS - ABORTERS - PRE STEP
    target = navigator.nextTarget(state, position, direction);
    BOOST_CHECK(!target.isNone());

    ACTS_INFO("<<< Test 1k-1o >>> step within 2nd layer at  "
              << toString(position));
  }

  // positive return: do the step
  step(position, direction, target);
  // (16) re-entering navigator:
  // POST STEP
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  // ACTORS - ABORTERS - PRE STEP
  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());

  ACTS_INFO("<<< Test 1p >>> step to 3rd layer at  " << toString(position));

  // Step through the surfaces on third layer
  for (std::size_t isf = 0; isf < 3; ++isf) {
    step(position, direction, target);
    // (17-19) re-entering navigator:
    // POST STEP
    navigator.handleSurfaceReached(state, position, direction,
                                   target.surface());
    // ACTORS - ABORTERS - PRE STEP
    target = navigator.nextTarget(state, position, direction);
    BOOST_CHECK(!target.isNone());

    ACTS_INFO("<<< Test 1q-1s >>> step within 3rd layer at  "
              << toString(position));
  }

  // positive return: do the step
  step(position, direction, target);
  // (20) re-entering navigator:
  // POST STEP
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  // ACTORS - ABORTERS - PRE STEP
  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());

  ACTS_INFO("<<< Test 1t >>> step to 4th layer at  " << toString(position));

  // Step through the surfaces on second layer
  for (std::size_t isf = 0; isf < 3; ++isf) {
    step(position, direction, target);
    // (21-23) re-entering navigator:
    // POST STEP
    navigator.handleSurfaceReached(state, position, direction,
                                   target.surface());
    // ACTORS - ABORTERS - PRE STEP
    target = navigator.nextTarget(state, position, direction);
    BOOST_CHECK(!target.isNone());

    ACTS_INFO("<<< Test 1t-1v >>> step within 4th layer at  "
              << toString(position));
  }

  // positive return: do the step
  step(position, direction, target);
  // (24) re-entering navigator:
  // POST STEP
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  // ACTORS - ABORTERS - PRE STEP
  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(target.isNone());

  ACTS_INFO("<<< Test 1w >>> step to boundary at  " << toString(position));
}

inline std::tuple<std::shared_ptr<const TrackingGeometry>,
                  std::vector<const Surface*>>
createDenseTelescope(const GeometryContext& geoCtx) {
  using namespace Acts;
  using namespace UnitLiterals;

  CuboidVolumeBuilder::Config conf;
  conf.position = {0., 0., 0.};
  conf.length = {2_m, 2_m, 2_m};

  {
    CuboidVolumeBuilder::SurfaceConfig surfaceTop;
    surfaceTop.position = {0, 0.5_m, 0.5_m};
    surfaceTop.rBounds = std::make_shared<RectangleBounds>(0.8_m, 0.2_m);

    CuboidVolumeBuilder::SurfaceConfig surfaceBottom;
    surfaceBottom.position = {0, -0.5_m, 0.5_m};
    surfaceBottom.rBounds = std::make_shared<RectangleBounds>(0.8_m, 0.2_m);

    CuboidVolumeBuilder::LayerConfig layer;
    layer.surfaceCfg.push_back(surfaceTop);
    layer.surfaceCfg.push_back(surfaceBottom);

    CuboidVolumeBuilder::VolumeConfig start;
    start.position = {0, 0, 0};
    start.length = {1.9_m, 1.9_m, 1.9_m};
    start.name = "start";
    start.layerCfg.push_back(layer);

    conf.volumeCfg.push_back(start);
  }

  CuboidVolumeBuilder cvb(conf);

  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });
  auto detector = TrackingGeometryBuilder(tgbCfg).trackingGeometry(geoCtx);

  std::vector<const Surface*> surfaces;
  detector->visitSurfaces(
      [&](const Surface* surface) { surfaces.push_back(surface); });

  return {std::move(detector), std::move(surfaces)};
}

BOOST_AUTO_TEST_CASE(Navigator_external_surfaces) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("NavigatorTest", logLevel));

  auto [detector, surfaces] = createDenseTelescope(tgContext);
  BOOST_CHECK_EQUAL(surfaces.size(), 2ul);
  const Surface& surfaceTop = *surfaces.at(0);
  const Surface& surfaceBottom = *surfaces.at(1);
  CHECK_CLOSE_ABS(surfaceTop.center(tgContext).y(), 0.5_m, 1e-6);
  CHECK_CLOSE_ABS(surfaceBottom.center(tgContext).y(), -0.5_m, 1e-6);

  Navigator::Config navCfg;
  navCfg.trackingGeometry = detector;
  navCfg.resolveSensitive = true;
  navCfg.resolveMaterial = true;
  navCfg.resolvePassive = false;
  Navigator navigator(navCfg, logger().clone("Navigator"));

  // check if we find no sensitive target starting from the middle without
  // external surfaces
  {
    ACTS_INFO("Test 1: start in the middle without external surfaces");

    Navigator::Options options(tgContext);
    Navigator::State state = navigator.makeState(options);

    Vector3 position = Vector3::Zero();
    Vector3 direction = Vector3::UnitZ();

    Result<void> result =
        navigator.initialize(state, position, direction, Direction::Forward());
    BOOST_CHECK(result.ok());

    NavigationTarget target = navigator.nextTarget(state, position, direction);

    BOOST_CHECK_NE(&target.surface(), surfaces.at(0));
    BOOST_CHECK_NE(&target.surface(), surfaces.at(1));
  }

  // check if we find a target starting from the top without external surfaces
  {
    ACTS_INFO("Test 2: start from top without external surfaces");

    Navigator::Options options(tgContext);
    Navigator::State state = navigator.makeState(options);

    Vector3 position = {0, 0.5_m, 0};
    Vector3 direction = Vector3::UnitZ();

    Result<void> result =
        navigator.initialize(state, position, direction, Direction::Forward());
    BOOST_CHECK(result.ok());

    NavigationTarget target = navigator.nextTarget(state, position, direction);

    BOOST_CHECK(!target.isNone());
    BOOST_CHECK_EQUAL(&target.surface(), &surfaceTop);
  }

  // check if we find a target starting from the bottom without external
  // surfaces
  {
    ACTS_INFO("Test 2: start from bottom without external surfaces");

    Navigator::Options options(tgContext);
    Navigator::State state = navigator.makeState(options);

    Vector3 position = {0, -0.5_m, 0};
    Vector3 direction = Vector3::UnitZ();

    Result<void> result =
        navigator.initialize(state, position, direction, Direction::Forward());
    BOOST_CHECK(result.ok());

    NavigationTarget target = navigator.nextTarget(state, position, direction);

    BOOST_CHECK(!target.isNone());
    BOOST_CHECK_EQUAL(&target.surface(), &surfaceBottom);
  }

  // check if we find the top surface starting from the middle with external
  // surfaces
  {
    ACTS_INFO("Test 3: start in the middle with external surfaces");

    Navigator::Options options(tgContext);
    options.insertExternalSurface(surfaceTop);
    Navigator::State state = navigator.makeState(options);

    Vector3 position = Vector3::Zero();
    Vector3 direction = Vector3::UnitZ();

    Result<void> result =
        navigator.initialize(state, position, direction, Direction::Forward());
    BOOST_CHECK(result.ok());

    NavigationTarget target = navigator.nextTarget(state, position, direction);

    BOOST_CHECK(!target.isNone());
    BOOST_CHECK_EQUAL(&target.surface(), &surfaceTop);
  }

  // check if we find the bottom surface starting from the top with external
  // surfaces
  {
    ACTS_INFO("Test 4: start from top with external surfaces");

    Navigator::Options options(tgContext);
    options.insertExternalSurface(surfaceBottom);
    Navigator::State state = navigator.makeState(options);

    Vector3 position = {0, 0.5_m, 0};
    Vector3 direction = Vector3::UnitZ();

    Result<void> result =
        navigator.initialize(state, position, direction, Direction::Forward());
    BOOST_CHECK(result.ok());

    NavigationTarget target = navigator.nextTarget(state, position, direction);

    BOOST_CHECK(!target.isNone());
    BOOST_CHECK_EQUAL(&target.surface(), &surfaceBottom);
  }

  // check if we find the top surface starting from the bottom with external
  // surfaces
  {
    ACTS_INFO("Test 5: start from bottom with external surfaces");

    Navigator::Options options(tgContext);
    options.insertExternalSurface(surfaceTop);
    Navigator::State state = navigator.makeState(options);

    Vector3 position = {0, -0.5_m, 0};
    Vector3 direction = Vector3::UnitZ();

    Result<void> result =
        navigator.initialize(state, position, direction, Direction::Forward());
    BOOST_CHECK(result.ok());

    NavigationTarget target = navigator.nextTarget(state, position, direction);

    BOOST_CHECK(!target.isNone());
    BOOST_CHECK_EQUAL(&target.surface(), &surfaceTop);
  }
}

BOOST_AUTO_TEST_CASE(TryAllNavigationPolicy_SurfaceInsideVolume) {
  auto logger = getDefaultLogger("UnitTests", Logging::VERBOSE);

  Experimental::Blueprint::Config cfg;
  cfg.envelope = ExtentEnvelope{{
      .z = {20_mm, 20_mm},
      .r = {0_mm, 20_mm},
  }};

  Experimental::Blueprint root{cfg};

  auto& cubcontainer =
      root.addCuboidContainer("CuboidContainer", AxisDirection::AxisZ);

  auto parentBounds = std::make_shared<CuboidVolumeBounds>(1_m, 1_m, 1_m);

  auto parentVol = std::make_unique<TrackingVolume>(Transform3::Identity(),
                                                    parentBounds, "parent");

  std::shared_ptr<const PlanarBounds> planarBounds =
      std::make_shared<const RectangleBounds>(5., 10.);

  auto surface = Surface::makeShared<PlaneSurface>(Transform3::Identity(),
                                                   std::move(planarBounds));

  auto detElement =
      std::make_unique<DetectorElementStub>(Transform3::Identity());

  surface->assignDetectorElement(*detElement);

  parentVol->assignGeometryId(GeometryIdentifier{}.withVolume(1));
  parentVol->addSurface(surface);
  auto parentNode =
      std::make_shared<Experimental::StaticBlueprintNode>(std::move(parentVol));

  // put two tracking volumes in the sides of the parent as children and a plane
  // surface in the middle of the parent volume
  double startZ1 = -1000. + 100. + 1.;
  double startZ2 = 1000. - 100. - 1.;
  Transform3 trf1 = Transform3(Translation3(0., 0., startZ1));
  Transform3 trf2 = Transform3(Translation3(0., 0., startZ2));

  auto childBounds = std::make_shared<CuboidVolumeBounds>(1_m, 1_m, 10_cm);
  auto childVol1 =
      std::make_unique<TrackingVolume>(trf1, childBounds, "child1");
  childVol1->assignGeometryId(GeometryIdentifier{}.withVolume(2));

  auto childNode1 =
      std::make_shared<Experimental::StaticBlueprintNode>(std::move(childVol1));

  auto childVol2 =
      std::make_unique<TrackingVolume>(trf2, childBounds, "child2");
  childVol2->assignGeometryId(GeometryIdentifier{}.withVolume(3));

  auto childNode2 =
      std::make_shared<Experimental::StaticBlueprintNode>(std::move(childVol2));

  parentNode->addChild(childNode1);
  parentNode->addChild(childNode2);

  cubcontainer.addChild(std::move(parentNode));

  auto trackingGeometry = root.construct({}, tgContext, *logger);

  Navigator::Config navCfg;
  navCfg.trackingGeometry =
      std::shared_ptr<const TrackingGeometry>(std::move(trackingGeometry));
  navCfg.resolveSensitive = true;
  navCfg.resolveMaterial = true;
  navCfg.resolvePassive = false;
  Navigator navigator(navCfg, logger->clone("Navigator"));

  Navigator::Options options(tgContext);
  Navigator::State state = navigator.makeState(options);

  Vector3 position{0., 0., -1000. + 0.5};
  Vector3 direction = Vector3::UnitZ();

  Result<void> result =
      navigator.initialize(state, position, direction, Direction::Forward());
  BOOST_CHECK(result.ok());
  BOOST_CHECK(state.currentVolume != nullptr);
  BOOST_CHECK(state.currentVolume->volumeName() == "parent");
  BOOST_CHECK_EQUAL(state.currentSurface, nullptr);

  // Do the chain simulating a straight line propagation
  // and check if we end up to the expected elements

  NavigationTarget target = navigator.nextTarget(state, position, direction);

  // it is supposed to find the boundary of the next child volume
  BOOST_CHECK(!target.isNone());
  auto targetGeoId = target.surface().geometryId();
  // portal of child1 volume as expected
  BOOST_CHECK(targetGeoId.volume() == 2 && targetGeoId.boundary() == 1);

  step(position, direction, target);
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  // check that we end up in the expected volume (child1)
  BOOST_CHECK(state.currentVolume->volumeName() == "child1");

  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());
  targetGeoId = target.surface().geometryId();
  // portal of child1 volume as expected
  BOOST_CHECK(targetGeoId.volume() == 2 && targetGeoId.boundary() == 2);

  step(position, direction, target);
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  // check that we end up in the expected volume (parent)
  BOOST_CHECK(state.currentVolume->volumeName() == "parent");
  // the next target should be the plane surface inside the parent volume
  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());
  BOOST_CHECK_EQUAL(&target.surface(), surface.get());

  step(position, direction, target);
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  // check that we are still in the parent volume
  BOOST_CHECK(state.currentVolume->volumeName() == "parent");

  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());
  targetGeoId = target.surface().geometryId();
  // portal of child2 volume as expected
  BOOST_CHECK(targetGeoId.volume() == 3 && targetGeoId.boundary() == 1);

  step(position, direction, target);
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  // check that we end up in the expected volume (child2)
  BOOST_CHECK(state.currentVolume->volumeName() == "child2");
  target = navigator.nextTarget(state, position, direction);
  BOOST_CHECK(!target.isNone());
  targetGeoId = target.surface().geometryId();
  // portal of child2 volume as expected
  BOOST_CHECK(targetGeoId.volume() == 3 && targetGeoId.boundary() == 2);

  step(position, direction, target);
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  // check that we end up in the expected volume (parent)
  BOOST_CHECK(state.currentVolume->volumeName() == "parent");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
