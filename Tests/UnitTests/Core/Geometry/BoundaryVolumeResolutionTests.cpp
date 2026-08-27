// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsTests/CommonHelpers/CubicTrackingGeometry.hpp"

#include <memory>
#include <optional>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

auto logger = getDefaultLogger("UnitTests", Logging::INFO);
auto gctx = GeometryContext::dangerouslyDefaultConstruct();

const TrackingVolume* findVolumeByName(const TrackingGeometry& trackingGeometry,
                                       const std::string& name) {
  const TrackingVolume* result = nullptr;
  trackingGeometry.apply([&](const TrackingVolume& volume) {
    if (volume.volumeName() == name) {
      result = &volume;
    }
  });
  return result;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GeometrySuite)

// Two cuboid volumes stacked in x with a shared (fused) portal at x=0. A
// position on the shared portal surface is ambiguous between the two
// volumes; the direction resolves the ambiguity.
BOOST_AUTO_TEST_CASE(BoundaryVolumeResolutionGen3) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisX] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisY] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  Blueprint root{cfg};

  root.addCuboidContainer("Stack", AxisDirection::AxisX, [&](auto& stack) {
    stack.addStaticVolume(
        Transform3{Translation3{Vector3{-100_mm, 0, 0}}},
        std::make_shared<CuboidVolumeBounds>(100_mm, 100_mm, 100_mm),
        "VolumeA");
    stack.addStaticVolume(
        Transform3{Translation3{Vector3{100_mm, 0, 0}}},
        std::make_shared<CuboidVolumeBounds>(100_mm, 100_mm, 100_mm),
        "VolumeB");
  });

  std::unique_ptr<TrackingGeometry> trackingGeometry =
      root.construct({}, gctx, *logger);
  BOOST_REQUIRE(trackingGeometry != nullptr);
  BOOST_CHECK(trackingGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen3);

  const TrackingVolume* volumeA =
      findVolumeByName(*trackingGeometry, "VolumeA");
  const TrackingVolume* volumeB =
      findVolumeByName(*trackingGeometry, "VolumeB");
  BOOST_REQUIRE(volumeA != nullptr);
  BOOST_REQUIRE(volumeB != nullptr);

  // The two volumes touch at x=0
  const Vector3 onBoundary{0, 0, 0};

  // Find the portal shared between the two volumes at the touching faces.
  // Note that the merged lateral portals of the stack are also shared
  // between the two volumes, so the position is needed to disambiguate.
  const Portal* sharedPortal = nullptr;
  for (const Portal& pa : volumeA->portals()) {
    for (const Portal& pb : volumeB->portals()) {
      if (&pa == &pb &&
          pa.surface().isOnSurface(gctx, onBoundary, Vector3::UnitX(),
                                   BoundaryTolerance::None())) {
        sharedPortal = &pa;
      }
    }
  }
  BOOST_REQUIRE(sharedPortal != nullptr);

  const Surface& portalSurface = sharedPortal->surface();

  // The direction-aware lookup resolves the correct side of the boundary
  auto forward = trackingGeometry->resolveLowestTrackingVolume(
      gctx, onBoundary, Vector3::UnitX(), &portalSurface);
  BOOST_REQUIRE(forward.ok());
  BOOST_CHECK_EQUAL(*forward, volumeB);

  auto backward = trackingGeometry->resolveLowestTrackingVolume(
      gctx, onBoundary, -Vector3::UnitX(), &portalSurface);
  BOOST_REQUIRE(backward.ok());
  BOOST_CHECK_EQUAL(*backward, volumeA);

  // Without direction and surface hint the lookup is position based and
  // returns one of the two adjacent volumes
  auto plain = trackingGeometry->resolveLowestTrackingVolume(gctx, onBoundary);
  BOOST_REQUIRE(plain.ok());
  BOOST_CHECK(*plain == volumeA || *plain == volumeB);
}

// Two glued Gen1 cuboid volumes sharing a boundary surface at x=0
BOOST_AUTO_TEST_CASE(BoundaryVolumeResolutionGen1) {
  CubicTrackingGeometry geometryBuilder{gctx};
  std::shared_ptr<const TrackingGeometry> trackingGeometry = geometryBuilder();
  BOOST_REQUIRE(trackingGeometry != nullptr);
  BOOST_CHECK(trackingGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen1);

  const TrackingVolume* volume1 =
      trackingGeometry->resolveLowestTrackingVolume(gctx, Vector3{-1.5_m, 0, 0})
          .value();
  const TrackingVolume* volume2 =
      trackingGeometry->resolveLowestTrackingVolume(gctx, Vector3{1.5_m, 0, 0})
          .value();
  BOOST_REQUIRE(volume1 != nullptr);
  BOOST_REQUIRE(volume2 != nullptr);
  BOOST_CHECK_EQUAL(volume1->volumeName(), "Volume 1");
  BOOST_CHECK_EQUAL(volume2->volumeName(), "Volume 2");

  // Find the glued boundary surface between the two volumes
  const Vector3 onBoundary{0, 0, 0};
  auto findBoundarySurface = [&](const TrackingVolume& volume,
                                 const Vector3& position) -> const Surface* {
    for (const auto& boundary : volume.boundarySurfaces()) {
      const Surface& surface = boundary->surfaceRepresentation();
      if (surface.isOnSurface(gctx, position, Vector3::UnitX(),
                              BoundaryTolerance::None())) {
        return &surface;
      }
    }
    return nullptr;
  };

  const Surface* boundarySurface = findBoundarySurface(*volume1, onBoundary);
  BOOST_REQUIRE(boundarySurface != nullptr);

  // The direction-aware lookup resolves the correct side of the boundary
  auto forward = trackingGeometry->resolveLowestTrackingVolume(
      gctx, onBoundary, Vector3::UnitX(), boundarySurface);
  BOOST_REQUIRE(forward.ok());
  BOOST_CHECK_EQUAL(*forward, volume2);

  auto backward = trackingGeometry->resolveLowestTrackingVolume(
      gctx, onBoundary, -Vector3::UnitX(), boundarySurface);
  BOOST_REQUIRE(backward.ok());
  BOOST_CHECK_EQUAL(*backward, volume1);

  // Outer boundary of volume 1 with nothing attached on the far side: this
  // is the end of the world
  const Vector3 onOuterBoundary{-3_m, 0, 0};
  const Surface* outerSurface = findBoundarySurface(*volume1, onOuterBoundary);
  BOOST_REQUIRE(outerSurface != nullptr);

  auto outward = trackingGeometry->resolveLowestTrackingVolume(
      gctx, onOuterBoundary, -Vector3::UnitX(), outerSurface);
  BOOST_REQUIRE(outward.ok());
  BOOST_CHECK_EQUAL(*outward, nullptr);

  // Pointing back inside stays in the volume
  auto inward = trackingGeometry->resolveLowestTrackingVolume(
      gctx, onOuterBoundary, Vector3::UnitX(), outerSurface);
  BOOST_REQUIRE(inward.ok());
  BOOST_CHECK_EQUAL(*inward, volume1);

  // A position slightly off the boundary is resolved through it within the
  // given tolerance
  const Vector3 nearBoundary{0.01_mm, 0, 0};
  auto nearResult = trackingGeometry->resolveLowestTrackingVolume(
      gctx, nearBoundary, -Vector3::UnitX(), boundarySurface, 0.1_mm);
  BOOST_REQUIRE(nearResult.ok());
  BOOST_CHECK_EQUAL(*nearResult, volume1);
}

// A position that is on the plane of a boundary surface of the resolved
// volume, but outside of its bounds, cannot be resolved through the boundary
BOOST_AUTO_TEST_CASE(BoundaryVolumeResolutionOffBounds) {
  CubicTrackingGeometry geometryBuilder{gctx};
  std::shared_ptr<const TrackingGeometry> trackingGeometry = geometryBuilder();
  BOOST_REQUIRE(trackingGeometry != nullptr);

  const TrackingVolume* volume1 =
      trackingGeometry->resolveLowestTrackingVolume(gctx, Vector3{-1.5_m, 0, 0})
          .value();
  BOOST_REQUIRE(volume1 != nullptr);

  // The glued boundary surface between the two volumes at x=0 spans
  // |y| < 0.5 m and |z| < 0.5 m
  const Surface* boundarySurface = nullptr;
  for (const auto& boundary : volume1->boundarySurfaces()) {
    const Surface& surface = boundary->surfaceRepresentation();
    if (surface.isOnSurface(gctx, Vector3::Zero(), Vector3::UnitX(),
                            BoundaryTolerance::None())) {
      boundarySurface = &surface;
    }
  }
  BOOST_REQUIRE(boundarySurface != nullptr);

  // Grazing the volume edge: within the lookup tolerance of the volume, but
  // outside of the boundary surface bounds
  const Vector3 offBounds{0, 0.5_m + 0.5_mm, 0};
  const double tolerance = 1_mm;
  BOOST_REQUIRE(trackingGeometry
                    ->resolveLowestTrackingVolume(gctx, offBounds, std::nullopt,
                                                  nullptr, tolerance)
                    .value() != nullptr);
  BOOST_REQUIRE(!boundarySurface->isOnSurface(
      gctx, offBounds, Vector3::UnitX(), BoundaryTolerance::None(), tolerance));

  auto result = trackingGeometry->resolveLowestTrackingVolume(
      gctx, offBounds, Vector3::UnitX(), boundarySurface, tolerance);
  BOOST_CHECK(!result.ok());
  BOOST_CHECK(result.error() == SurfaceError::GlobalPositionNotOnSurface);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
