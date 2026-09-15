// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"

#include <cmath>
#include <stdexcept>

#include "NavigatorTelescope.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace ActsTests::NavigatorTelescope;

namespace ActsTests {

// `NavigatorPlainOptions::boundaryToleranceOverrides` relaxes the bounds check
// on a surface of the tracking geometry. The geometry still resolves the
// surface, in its volume for Gen3 and on its layer for Gen1, so both
// generations are covered.

BOOST_AUTO_TEST_SUITE(NavigatorBoundaryToleranceOverride)

// Without an override the navigator must not target a surface the track misses

BOOST_AUTO_TEST_CASE(BaselineGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  BOOST_REQUIRE(telescope.geometry->geometryVersion() ==
                TrackingGeometry::GeometryVersion::Gen1);
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  Navigator::Options options(gctx);
  NavigationTarget target = firstTarget(navigator, options, Vector3::Zero());

  BOOST_CHECK_NE(&target.surface(), telescope.top);
  BOOST_CHECK_NE(&target.surface(), telescope.bottom);
}

BOOST_AUTO_TEST_CASE(BaselineGen3) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  BOOST_REQUIRE(telescope.geometry->geometryVersion() ==
                TrackingGeometry::GeometryVersion::Gen3);
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  Navigator::Options options(gctx);
  NavigationTarget target = firstTarget(navigator, options, Vector3::Zero());

  BOOST_CHECK_NE(&target.surface(), telescope.top);
  BOOST_CHECK_NE(&target.surface(), telescope.bottom);
}

// The bounds check is dropped, so the surface the track misses is targeted

BOOST_AUTO_TEST_CASE(BoundsExtensionGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  Navigator::Options options(gctx);
  options.overrideBoundaryTolerance(*telescope.top);
  NavigationTarget target = firstTarget(navigator, options, Vector3::Zero());

  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK_EQUAL(&target.surface(), telescope.top);
}

// The policy returns the same surface with a `None` tolerance. The override
// goes in first, so the de-duplication keeps the tolerance of the caller.
BOOST_AUTO_TEST_CASE(BoundsExtensionGen3) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  Navigator::Options options(gctx);
  options.overrideBoundaryTolerance(*telescope.top);
  NavigationTarget target = firstTarget(navigator, options, Vector3::Zero());

  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK_EQUAL(&target.surface(), telescope.top);
  BOOST_CHECK(!target.isPortalTarget());
}

// The tolerance is per entry, so a caller can ask for a bounds check

BOOST_AUTO_TEST_CASE(BoundsCheckedGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  Navigator::Options options(gctx);
  options.overrideBoundaryTolerance(*telescope.top, BoundaryTolerance::None());
  NavigationTarget target = firstTarget(navigator, options, Vector3::Zero());

  BOOST_CHECK_NE(&target.surface(), telescope.top);
}

BOOST_AUTO_TEST_CASE(BoundsCheckedGen3) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  Navigator::Options options(gctx);
  options.overrideBoundaryTolerance(*telescope.top, BoundaryTolerance::None());
  NavigationTarget target = firstTarget(navigator, options, Vector3::Zero());

  BOOST_CHECK_NE(&target.surface(), telescope.top);
}

// The geometry has to hold the surface, otherwise it never resolves it and the
// override is silently lost. The check is shared by both generations.

BOOST_AUTO_TEST_CASE(SurfaceOutsideTheGeometryThrows) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto outOfGeometry = makeOutOfGeometrySurface({0, 0, 0.25_m});
  BOOST_REQUIRE(outOfGeometry->geometryId() == GeometryIdentifier{});

  Navigator::Options options(gctx);
  options.overrideBoundaryTolerance(*outOfGeometry);
  Navigator::State state = navigator.makeState(options);
  Vector3 position = Vector3::Zero();
  const Vector3 direction = Vector3::UnitZ();
  BOOST_CHECK_THROW(static_cast<void>(navigator.initialize(
                        state, position, direction, Direction::Forward())),
                    std::invalid_argument);
}

// The relaxed bounds check must not pull the propagation out of its volume.
// The track below crosses the plane of the top surface at y = 1.2 m, outside
// the 0.95 m volume, and leaves through the boundary at y = 0.95 m first.

namespace {
const Vector3 pokeStart{0, 0.9_m, 0};
const Vector3 pokeDir = Vector3{0, 0.6, 1}.normalized();
}  // namespace

// Gen1 resolves the surfaces of a layer before the boundaries, so it checks
// the intersection against the volume explicitly.
BOOST_AUTO_TEST_CASE(CrossingOutsideTheVolumeIsDroppedGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  Navigator::Options options(gctx);
  options.overrideBoundaryTolerance(*telescope.top);
  Navigator::State state = navigator.makeState(options);
  Vector3 position = pokeStart;
  BOOST_REQUIRE(
      navigator.initialize(state, position, pokeDir, Direction::Forward())
          .ok());

  NavigationTarget target = navigator.nextTarget(state, position, pokeDir);
  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK(target.surface().geometryId().boundary() != 0);
}

// Gen3 sorts portals and surfaces together, so the portal wins on its own.
BOOST_AUTO_TEST_CASE(CrossingOutsideTheVolumeIsDroppedGen3) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  Navigator::Options options(gctx);
  options.overrideBoundaryTolerance(*telescope.top);
  Navigator::State state = navigator.makeState(options);
  Vector3 position = pokeStart;
  BOOST_REQUIRE(
      navigator.initialize(state, position, pokeDir, Direction::Forward())
          .ok());

  NavigationTarget target = navigator.nextTarget(state, position, pokeDir);
  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK(target.isPortalTarget());
}

// The override applies where the geometry resolves the surface: in its volume
// for Gen3, on its layer for Gen1

BOOST_AUTO_TEST_CASE(ScopedToItsVolumeGen3) {
  auto logger = getDefaultLogger("Gen3", logLevel);

  Blueprint::Config cfg;
  cfg.envelope = ExtentEnvelope{{
      .x = {20_mm, 20_mm},
      .y = {20_mm, 20_mm},
      .z = {20_mm, 20_mm},
  }};
  Blueprint root{cfg};

  // Two volumes along z. The far one holds a small surface off the axis.
  const Surface* overridden = nullptr;
  root.addCuboidContainer("Stack", AxisDirection::AxisZ, [&](auto& stack) {
    stack.addChild(
        std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
            Transform3{Translation3{Vector3{0, 0, -0.5_m}}},
            std::make_shared<CuboidVolumeBounds>(0.5_m, 0.5_m, 0.5_m),
            "near")));

    auto farVolume = std::make_unique<TrackingVolume>(
        Transform3{Translation3{Vector3{0, 0, 0.5_m}}},
        std::make_shared<CuboidVolumeBounds>(0.5_m, 0.5_m, 0.5_m), "far");
    auto surface = Surface::makeShared<PlaneSurface>(
        Transform3{Translation3{Vector3{0, 0.4_m, 0.5_m}}},
        std::make_shared<const RectangleBounds>(0.05_m, 0.05_m));
    overridden = surface.get();
    farVolume->addSurface(std::move(surface));
    stack.addChild(std::make_shared<StaticBlueprintNode>(std::move(farVolume)));
  });

  auto geometry = root.construct({}, gctx, *logger);
  BOOST_REQUIRE(overridden != nullptr);
  Navigator navigator = makeNavigator(
      std::shared_ptr<const TrackingGeometry>(std::move(geometry)), *logger);

  Navigator::Options options(gctx);
  options.overrideBoundaryTolerance(*overridden);
  Navigator::State state = navigator.makeState(options);

  Vector3 position{0, 0, -0.9_m};
  const Vector3 dir = Vector3::UnitZ();
  BOOST_REQUIRE(
      navigator.initialize(state, position, dir, Direction::Forward()).ok());

  // The navigator resolved the volume that holds the surface
  BOOST_REQUIRE_EQUAL(state.boundaryToleranceOverrides.size(), 1u);
  BOOST_REQUIRE(state.boundaryToleranceOverrides.front().volume != nullptr);
  BOOST_CHECK_EQUAL(
      state.boundaryToleranceOverrides.front().volume->volumeName(), "far");

  // In the near volume the portal is targeted, not the overridden surface
  NavigationTarget target = navigator.nextTarget(state, position, dir);
  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK(target.isPortalTarget());
  stepOnto(position, dir, target.surface());
  navigator.handleSurfaceReached(state, position, dir, target.surface());
  BOOST_REQUIRE_EQUAL(state.currentVolume->volumeName(), "far");

  // In its own volume the relaxed bounds check applies and it is targeted
  target = navigator.nextTarget(state, position, dir);
  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK_EQUAL(&target.surface(), overridden);
}

BOOST_AUTO_TEST_CASE(ScopedToItsLayerGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));

  CuboidVolumeBuilder::Config conf;
  conf.position = {0., 0., 0.};
  conf.length = {2_m, 2_m, 2_m};

  // One layer on the axis, a second one further out, off the axis
  CuboidVolumeBuilder::SurfaceConfig onAxis;
  onAxis.position = {0, 0, -0.3_m};
  onAxis.rBounds = std::make_shared<RectangleBounds>(0.1_m, 0.1_m);
  CuboidVolumeBuilder::LayerConfig nearLayer;
  nearLayer.binningDimension = AxisDirection::AxisZ;
  nearLayer.surfaceCfg.push_back(onAxis);

  CuboidVolumeBuilder::SurfaceConfig offAxis;
  offAxis.position = {0, 0.4_m, 0.3_m};
  offAxis.rBounds = std::make_shared<RectangleBounds>(0.05_m, 0.05_m);
  CuboidVolumeBuilder::LayerConfig farLayer;
  farLayer.binningDimension = AxisDirection::AxisZ;
  farLayer.surfaceCfg.push_back(offAxis);

  CuboidVolumeBuilder::VolumeConfig volume;
  volume.binningDimension = AxisDirection::AxisZ;
  volume.position = {0, 0, 0};
  volume.length = {1.9_m, 1.9_m, 1.9_m};
  volume.name = "telescope";
  volume.layerCfg.push_back(nearLayer);
  volume.layerCfg.push_back(farLayer);
  conf.volumeCfg.push_back(volume);

  CuboidVolumeBuilder cvb(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });
  auto geometry = TrackingGeometryBuilder(tgbCfg).trackingGeometry(gctx);

  std::vector<const Surface*> surfaces;
  geometry->visitSurfaces([&](const Surface* sf) { surfaces.push_back(sf); });
  BOOST_REQUIRE_EQUAL(surfaces.size(), 2u);
  const Surface* overridden = surfaces.at(1);
  BOOST_REQUIRE_EQUAL(overridden->center(gctx).y(), 0.4_m);

  Navigator navigator = makeNavigator(
      std::shared_ptr<const TrackingGeometry>(std::move(geometry)), logger());

  Navigator::Options options(gctx);
  options.overrideBoundaryTolerance(*overridden);
  std::vector<const Surface*> reached =
      walk(navigator, options, Vector3{0, 0, -0.8_m}, Vector3::UnitZ());

  // The override applies only while its layer is the current one, so the near
  // layer comes first
  BOOST_REQUIRE(!reached.empty());
  BOOST_CHECK_NE(reached.front(), overridden);
  BOOST_CHECK_LT(
      std::abs(reached.front()->center(gctx).z() - onAxis.position.z()), 1_um);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
