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
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

#include <algorithm>

#include "NavigatorTelescope.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace ActsTests::NavigatorTelescope;

namespace ActsTests {

// Tests of `NavigatorPlainOptions::additionalSurfaces` for Gen1 and Gen3

namespace {

/// Whether @p sub appears in @p full in the same order, gaps allowed.
bool isSubsequence(const std::vector<const Surface*>& sub,
                   const std::vector<const Surface*>& full) {
  auto itr = full.begin();
  for (const Surface* surface : sub) {
    itr = std::find(itr, full.end(), surface);
    if (itr == full.end()) {
      return false;
    }
    ++itr;
  }
  return true;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(NavigatorAdditionalSurface)

// By default a surface the track misses is targeted too

BOOST_AUTO_TEST_CASE(OffPathGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto additional = makeOutOfGeometrySurface({0, surfaceY, 0.25_m});

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*additional);
  NavigationTarget target = firstTarget(navigator, options, Vector3::Zero());

  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK_EQUAL(&target.surface(), additional.get());
}

BOOST_AUTO_TEST_CASE(OffPathGen3) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto additional = makeOutOfGeometrySurface({0, surfaceY, 0.25_m});

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*additional);
  NavigationTarget target = firstTarget(navigator, options, Vector3::Zero());

  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK_EQUAL(&target.surface(), additional.get());
}

BOOST_AUTO_TEST_CASE(BoundsCheckedGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto additional = makeOutOfGeometrySurface({0, surfaceY, 0.25_m});

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*additional, BoundaryTolerance::None());
  std::vector<const Surface*> reached =
      walk(navigator, options, Vector3::Zero(), Vector3::UnitZ());

  BOOST_CHECK_EQUAL(std::ranges::count(reached, additional.get()), 0);
}

BOOST_AUTO_TEST_CASE(BoundsCheckedGen3) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto additional = makeOutOfGeometrySurface({0, surfaceY, 0.25_m});

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*additional, BoundaryTolerance::None());
  std::vector<const Surface*> reached =
      walk(navigator, options, Vector3::Zero(), Vector3::UnitZ());

  BOOST_CHECK_EQUAL(std::ranges::count(reached, additional.get()), 0);
}

// The candidates of the geometry are all still reached

BOOST_AUTO_TEST_CASE(NothingOfTheGeometryIsSkippedGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  // Between the origin and the two telescope surfaces
  auto additional = makeOutOfGeometrySurface({0, 0, 0.25_m});

  Navigator::Options plain(gctx);
  std::vector<const Surface*> baseline =
      walk(navigator, plain, Vector3::Zero(), Vector3::UnitZ(), 20);

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*additional);
  std::vector<const Surface*> reached =
      walk(navigator, options, Vector3::Zero(), Vector3::UnitZ(), 20);

  BOOST_REQUIRE(!baseline.empty());
  BOOST_CHECK_EQUAL(std::ranges::count(reached, additional.get()), 1);
  BOOST_CHECK(isSubsequence(baseline, reached));
}

BOOST_AUTO_TEST_CASE(NothingOfTheGeometryIsSkippedGen3) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto additional = makeOutOfGeometrySurface({0, 0, 0.25_m});

  Navigator::Options plain(gctx);
  std::vector<const Surface*> baseline =
      walk(navigator, plain, Vector3::Zero(), Vector3::UnitZ(), 20);

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*additional);
  std::vector<const Surface*> reached =
      walk(navigator, options, Vector3::Zero(), Vector3::UnitZ(), 20);

  BOOST_REQUIRE(!baseline.empty());
  BOOST_CHECK_EQUAL(std::ranges::count(reached, additional.get()), 1);
  BOOST_CHECK(isSubsequence(baseline, reached));
}

// The surface crosses the track at z = 0.96 m, outside the telescope volume,
// so the closer boundary wins.
BOOST_AUTO_TEST_CASE(PortalWinsGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto outside = makeOutOfGeometrySurface({0, 0, 0.96_m}, 10_m, 10_m);

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*outside);
  std::vector<const Surface*> reached =
      walk(navigator, options, Vector3::Zero(), Vector3::UnitZ(), 20);

  auto boundary = std::ranges::find_if(reached, [](const Surface* sf) {
    return sf->geometryId().boundary() != 0;
  });
  BOOST_REQUIRE(boundary != reached.end());
  BOOST_CHECK_EQUAL(std::count(reached.begin(), boundary, outside.get()), 0);
}

// The Gen3 world volume reaches z = 0.97 m, so the surface is reached there.
BOOST_AUTO_TEST_CASE(PortalWinsGen3) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto outside = makeOutOfGeometrySurface({0, 0, 0.96_m}, 10_m, 10_m);

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*outside);
  std::vector<const Surface*> reached =
      walk(navigator, options, Vector3::Zero(), Vector3::UnitZ(), 20);

  auto itr = std::ranges::find(reached, outside.get());
  BOOST_REQUIRE(itr != reached.end());
  BOOST_CHECK_EQUAL(std::ranges::count(reached, outside.get()), 1);
  // The portal of the telescope volume comes first
  BOOST_CHECK(std::ranges::any_of(reached.begin(), itr, [](const Surface* sf) {
    return sf->geometryId().boundary() != 0;
  }));
}

// A line surface, intersected at the point of closest approach
BOOST_AUTO_TEST_CASE(PerigeeGen3) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto perigee = Surface::makeShared<PerigeeSurface>(Vector3{0.1_m, 0, 0.25_m});

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*perigee);
  const Vector3 dir = Vector3{0.3, 0, 1}.normalized();
  std::vector<const Surface*> reached =
      walk(navigator, options, Vector3::Zero(), dir);

  BOOST_CHECK_EQUAL(std::ranges::count(reached, perigee.get()), 1);
}

// `dropAfterReached`. A straight walk cannot approach a surface twice, so the
// propagation steps back in front of it.

namespace {
/// Reach @p additional, then ask for the next target from just before it.
NavigationTarget targetOnSecondApproach(const Navigator& navigator,
                                        const Navigator::Options& options,
                                        const Surface& additional) {
  Navigator::State state = navigator.makeState(options);
  Vector3 position = Vector3::Zero();
  const Vector3 direction = Vector3::UnitZ();
  BOOST_REQUIRE(
      navigator
          .initialize(state, {.position = position, .direction = direction})
          .ok());

  NavigationTarget target = navigator.nextTarget(state, position, direction);
  BOOST_REQUIRE(!target.isNone());
  BOOST_REQUIRE_EQUAL(&target.surface(), &additional);
  stepOnto(position, direction, target.surface());
  navigator.handleSurfaceReached(state, position, direction, target.surface());

  // Approach it once more, from 1 mm in front of it
  position -= 1_mm * direction;
  return navigator.nextTarget(state, position, direction);
}
}  // namespace

BOOST_AUTO_TEST_CASE(DroppedAfterReachedGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto additional = makeOutOfGeometrySurface({0, 0, 0.25_m});

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*additional);
  NavigationTarget target =
      targetOnSecondApproach(navigator, options, *additional);

  BOOST_CHECK_NE(&target.surface(), additional.get());
}

BOOST_AUTO_TEST_CASE(KeptAfterReachedGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto additional = makeOutOfGeometrySurface({0, 0, 0.25_m});

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*additional, BoundaryTolerance::Infinite(),
                                    nullptr, false);
  NavigationTarget target =
      targetOnSecondApproach(navigator, options, *additional);

  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK_EQUAL(&target.surface(), additional.get());
}

// A volume scopes the surface to that one volume

BOOST_AUTO_TEST_CASE(ScopedToAVolumeGen3) {
  auto logger = getDefaultLogger("Gen3", logLevel);

  Blueprint::Config cfg;
  cfg.envelope = ExtentEnvelope{{
      .x = {20_mm, 20_mm},
      .y = {20_mm, 20_mm},
      .z = {20_mm, 20_mm},
  }};
  Blueprint root{cfg};

  root.addCuboidContainer("Stack", AxisDirection::AxisZ, [&](auto& stack) {
    stack.addChild(
        std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
            Transform3{Translation3{Vector3{0, 0, -0.5_m}}},
            std::make_shared<CuboidVolumeBounds>(0.5_m, 0.5_m, 0.5_m),
            "near")));
    stack.addChild(
        std::make_shared<StaticBlueprintNode>(std::make_unique<TrackingVolume>(
            Transform3{Translation3{Vector3{0, 0, 0.5_m}}},
            std::make_shared<CuboidVolumeBounds>(0.5_m, 0.5_m, 0.5_m), "far")));
  });

  auto geometry = std::shared_ptr<const TrackingGeometry>(
      root.construct({}, gctx, *logger));
  Navigator navigator = makeNavigator(geometry, *logger);

  // Inside the near volume
  auto additional = makeOutOfGeometrySurface({0, 0, -0.5_m});
  const TrackingVolume* farVolume = geometry->findVolumeByName("far");
  BOOST_REQUIRE(farVolume != nullptr);

  const Vector3 start{0, 0, -0.9_m};

  // Scoped to the far volume, where the track never crosses it
  {
    Navigator::Options options(gctx);
    options.registerAdditionalSurface(*additional,
                                      BoundaryTolerance::Infinite(), farVolume);
    std::vector<const Surface*> reached =
        walk(navigator, options, start, Vector3::UnitZ(), 12);
    BOOST_CHECK_EQUAL(std::ranges::count(reached, additional.get()), 0);
  }
  // Unscoped, so the volume that holds the crossing offers it
  {
    Navigator::Options options(gctx);
    options.registerAdditionalSurface(*additional);
    std::vector<const Surface*> reached =
        walk(navigator, options, start, Vector3::UnitZ(), 12);
    BOOST_CHECK_EQUAL(std::ranges::count(reached, additional.get()), 1);
  }
}

// Several additional surfaces are offered in the order the track crosses them

BOOST_AUTO_TEST_CASE(SeveralSurfacesInOrderGen3) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen3", logLevel));
  Telescope telescope = makeTelescopeGen3(logger());
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  auto first = makeOutOfGeometrySurface({0, 0, 0.1_m});
  auto second = makeOutOfGeometrySurface({0, 0, 0.2_m});
  auto third = makeOutOfGeometrySurface({0, 0, 0.3_m});

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*third);
  options.registerAdditionalSurface(*first);
  options.registerAdditionalSurface(*second);

  std::vector<const Surface*> reached =
      walk(navigator, options, Vector3::Zero(), Vector3::UnitZ(), 20);

  std::vector<const Surface*> reachedAdditional;
  for (const Surface* surface : reached) {
    if (surface == first.get() || surface == second.get() ||
        surface == third.get()) {
      reachedAdditional.push_back(surface);
    }
  }
  const std::vector<const Surface*> expected{first.get(), second.get(),
                                             third.get()};
  BOOST_CHECK(reachedAdditional == expected);
}

// The closer solution is behind the propagation, so the one ahead is offered

BOOST_AUTO_TEST_CASE(SecondSolutionAheadGen1) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("Gen1", logLevel));
  Telescope telescope = makeTelescopeGen1();
  Navigator navigator = makeNavigator(telescope.geometry, logger());

  // A cylinder around the y axis, so a track along z crosses it at z = +-R
  constexpr double radius = 0.3_m;
  auto additional = Surface::makeShared<CylinderSurface>(
      Transform3{AngleAxis3{90._degree, Vector3::UnitX()}},
      std::make_shared<const CylinderBounds>(radius, 0.5_m));

  Navigator::Options options(gctx);
  options.registerAdditionalSurface(*additional);
  // Start closer to the solution behind than to the one ahead
  NavigationTarget target =
      firstTarget(navigator, options, {0, 0, -0.2 * radius});

  BOOST_REQUIRE(!target.isNone());
  BOOST_CHECK_EQUAL(&target.surface(), additional.get());
  BOOST_CHECK_GT(target.pathLength(), 0.);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
