// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <vector>

/// A telescope of two plane surfaces, built for both geometry generations.
/// Shared by the navigator tests of the boundary tolerance overrides and of
/// the external surfaces.
namespace ActsTests::NavigatorTelescope {

using namespace Acts;
using namespace Acts::UnitLiterals;

const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

const Logging::Level logLevel = Logging::INFO;

// One 2x2x2 m volume around the origin with two plane surfaces at z = 0.5 m,
// displaced to y = +-0.5 m. Their half extent in y is 0.2 m, so a track along
// z through the origin misses both by 0.3 m.
constexpr double surfaceZ = 0.5_m;
constexpr double surfaceY = 0.5_m;
constexpr double surfaceHalfX = 0.8_m;
constexpr double surfaceHalfY = 0.2_m;

struct Telescope {
  std::shared_ptr<const TrackingGeometry> geometry;
  const Surface* top{};
  const Surface* bottom{};
};

/// Build the telescope as a Gen1 tracking geometry.
inline Telescope makeTelescopeGen1() {
  CuboidVolumeBuilder::Config conf;
  conf.position = {0., 0., 0.};
  conf.length = {2_m, 2_m, 2_m};

  CuboidVolumeBuilder::SurfaceConfig top;
  top.position = {0, surfaceY, surfaceZ};
  top.rBounds = std::make_shared<RectangleBounds>(surfaceHalfX, surfaceHalfY);

  CuboidVolumeBuilder::SurfaceConfig bottom;
  bottom.position = {0, -surfaceY, surfaceZ};
  bottom.rBounds =
      std::make_shared<RectangleBounds>(surfaceHalfX, surfaceHalfY);

  CuboidVolumeBuilder::LayerConfig layer;
  layer.binningDimension = AxisDirection::AxisZ;
  layer.surfaceCfg.push_back(top);
  layer.surfaceCfg.push_back(bottom);

  CuboidVolumeBuilder::VolumeConfig volume;
  volume.binningDimension = AxisDirection::AxisZ;
  volume.position = {0, 0, 0};
  volume.length = {1.9_m, 1.9_m, 1.9_m};
  volume.name = "telescope";
  volume.layerCfg.push_back(layer);

  conf.volumeCfg.push_back(volume);

  CuboidVolumeBuilder cvb(conf);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });

  Telescope telescope;
  telescope.geometry = TrackingGeometryBuilder(tgbCfg).trackingGeometry(gctx);

  std::vector<const Surface*> surfaces;
  telescope.geometry->visitSurfaces(
      [&](const Surface* surface) { surfaces.push_back(surface); });
  BOOST_REQUIRE_EQUAL(surfaces.size(), 2u);
  telescope.top = surfaces.at(0);
  telescope.bottom = surfaces.at(1);
  return telescope;
}

/// Build the same telescope as a Gen3 tracking geometry.
inline Telescope makeTelescopeGen3(const Logger& logger) {
  Blueprint::Config cfg;
  cfg.envelope = ExtentEnvelope{{
      .x = {20_mm, 20_mm},
      .y = {20_mm, 20_mm},
      .z = {20_mm, 20_mm},
  }};
  Blueprint root{cfg};

  auto volume = std::make_unique<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(0.95_m, 0.95_m, 0.95_m),
      "telescope");

  Telescope telescope;
  for (double y : {surfaceY, -surfaceY}) {
    auto surface = Surface::makeShared<PlaneSurface>(
        Transform3{Translation3{Vector3{0, y, surfaceZ}}},
        std::make_shared<const RectangleBounds>(surfaceHalfX, surfaceHalfY));
    if (y > 0) {
      telescope.top = surface.get();
    } else {
      telescope.bottom = surface.get();
    }
    volume->addSurface(std::move(surface));
  }

  root.addChild(std::make_shared<StaticBlueprintNode>(std::move(volume)));
  telescope.geometry = root.construct({}, gctx, logger);
  return telescope;
}

/// A plane surface with its normal along z, held by no volume, so it keeps
/// the default geometry identifier.
inline std::shared_ptr<Surface> makeOutOfGeometrySurface(const Vector3& center,
                                                         double halfX = 0.1_m,
                                                         double halfY = 0.1_m) {
  return Surface::makeShared<PlaneSurface>(
      Transform3{Translation3{center}},
      std::make_shared<const RectangleBounds>(halfX, halfY));
}

/// Move @p position onto @p surface along @p direction
inline void stepOnto(Vector3& position, const Vector3& direction,
                     const Surface& surface) {
  const Intersection3D intersection =
      surface.intersect(gctx, position, direction).closestForward();
  BOOST_REQUIRE(intersection.isValid());
  position += intersection.pathLength() * direction;
}

/// Initialize the navigator at @p position along +z, return its first target
inline NavigationTarget firstTarget(const Navigator& navigator,
                                    const Navigator::Options& options,
                                    const Vector3& position) {
  Navigator::State state = navigator.makeState(options);
  Vector3 pos = position;
  const Vector3 dir = Vector3::UnitZ();
  Result<void> result =
      navigator.initialize(state, pos, dir, Direction::Forward());
  BOOST_REQUIRE(result.ok());
  return navigator.nextTarget(state, pos, dir);
}

/// Walk the navigation and return every surface the propagation reaches.
inline std::vector<const Surface*> walk(const Navigator& navigator,
                                        const Navigator::Options& options,
                                        const Vector3& start,
                                        const Vector3& dir, int maxSteps = 12) {
  Navigator::State state = navigator.makeState(options);
  Vector3 position = start;
  BOOST_REQUIRE(
      navigator.initialize(state, position, dir, Direction::Forward()).ok());

  std::vector<const Surface*> reached;
  for (int i = 0; i < maxSteps; ++i) {
    NavigationTarget target = navigator.nextTarget(state, position, dir);
    if (target.isNone()) {
      break;
    }
    reached.push_back(&target.surface());
    stepOnto(position, dir, target.surface());
    navigator.handleSurfaceReached(state, position, dir, target.surface());
  }
  return reached;
}

inline Navigator makeNavigator(std::shared_ptr<const TrackingGeometry> geometry,
                               const Logger& logger) {
  Navigator::Config cfg;
  cfg.trackingGeometry = std::move(geometry);
  cfg.resolveSensitive = true;
  cfg.resolveMaterial = true;
  cfg.resolvePassive = true;
  return Navigator(cfg, logger.clone("Navigator"));
}

}  // namespace ActsTests::NavigatorTelescope
