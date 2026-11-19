// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumePlacementBase.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

const GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

constexpr double rMin = 10_mm;
constexpr double rMax = 20_mm;
constexpr double halfZ = 30_mm;

/// Placement with a fixed transform, which makes a volume alignable
class FixedPlacement final : public VolumePlacementBase {
 public:
  explicit FixedPlacement(const Transform3& transform)
      : m_transform{transform}, m_inverse{transform.inverse()} {}

  const Transform3& localToGlobalTransform(
      const GeometryContext& /*gctx*/) const override {
    return m_transform;
  }

  const Transform3& globalToLocalTransform(
      const GeometryContext& /*gctx*/) const override {
    return m_inverse;
  }

  const Transform3& portalLocalToGlobal(const GeometryContext& /*gctx*/,
                                        std::size_t portalIdx) const override {
    return m_portals.at(portalIdx);
  }

  void makePortalsAlignable(const GeometryContext& context,
                            const std::vector<std::shared_ptr<RegularSurface>>&
                                portalsToAlign) override {
    VolumePlacementBase::makePortalsAlignable(context, portalsToAlign);
    for (std::size_t i = 0; i < portalsToAlign.size(); ++i) {
      m_portals.push_back(alignPortal(context, i));
    }
  }

 private:
  Transform3 m_transform;
  Transform3 m_inverse;
  std::vector<Transform3> m_portals;
};

/// Three cylinders stacked along z, so their inner and outer cylinder portals
/// are merged into portals three times as long as each volume.
std::shared_ptr<TrackingGeometry> makeCylinderStack(const Logger& logger) {
  Blueprint::Config cfg;
  cfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  cfg.envelope[AxisDirection::AxisR] = {2_mm, 20_mm};
  Blueprint root{cfg};

  auto& stack = root.addCylinderContainer("Stack", AxisDirection::AxisZ);
  auto bounds = std::make_shared<CylinderVolumeBounds>(rMin, rMax, halfZ);
  for (int i = 0; i < 3; ++i) {
    stack.addStaticVolume(std::make_unique<TrackingVolume>(
        Transform3{Translation3{Vector3{0, 0, (i - 1) * 2 * halfZ}}}, bounds,
        "child" + std::to_string(i)));
  }
  return root.construct({}, gctx, logger);
}

const TrackingVolume& findVolume(const TrackingGeometry& geometry,
                                 std::string_view name) {
  const TrackingVolume* found = nullptr;
  geometry.visitVolumes([&](const TrackingVolume* volume) {
    if (volume->volumeName() == name) {
      found = volume;
    }
  });
  BOOST_REQUIRE(found != nullptr);
  return *found;
}

const Portal& outerCylinder(const TrackingVolume& volume) {
  for (const Portal& portal : volume.portals()) {
    const auto* bounds =
        dynamic_cast<const CylinderBounds*>(&portal.surface().bounds());
    if (bounds != nullptr &&
        std::abs(bounds->get(CylinderBounds::eR) - rMax) < 1e-6) {
      return portal;
    }
  }
  BOOST_FAIL("no outer cylinder portal");
  throw std::logic_error("unreachable");
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagatorSuite)
BOOST_AUTO_TEST_SUITE(NavigatorPortalFace)

// Merged portals know the face of each volume, others need none

BOOST_AUTO_TEST_CASE(CylinderStackFaces) {
  auto logger = getDefaultLogger("PortalFace", Logging::INFO);
  auto geometry = makeCylinderStack(*logger);

  for (int i = 0; i < 3; ++i) {
    const TrackingVolume& volume =
        findVolume(*geometry, "child" + std::to_string(i));
    for (const Portal& portal : volume.portals()) {
      const Surface* face = volume.portalFace(portal);
      if (portal.surface().type() != Surface::Cylinder) {
        // Discs between the volumes are fused, the outer ones are not merged
        BOOST_CHECK(face == nullptr);
        continue;
      }
      BOOST_REQUIRE(face != nullptr);
      BOOST_CHECK(face != &portal.surface());
      BOOST_CHECK_EQUAL(face->type(), Surface::Cylinder);
      const auto& faceBounds =
          dynamic_cast<const CylinderBounds&>(face->bounds());
      const auto& portalBounds =
          dynamic_cast<const CylinderBounds&>(portal.surface().bounds());
      BOOST_CHECK_CLOSE(faceBounds.get(CylinderBounds::eR),
                        portalBounds.get(CylinderBounds::eR), 1e-9);
      BOOST_CHECK_CLOSE(faceBounds.get(CylinderBounds::eHalfLengthZ), halfZ,
                        1e-9);
      BOOST_CHECK_CLOSE(portalBounds.get(CylinderBounds::eHalfLengthZ),
                        3 * halfZ, 1e-9);
      BOOST_CHECK_SMALL((face->center(gctx) - volume.center(gctx)).norm(),
                        1e-9);
    }
  }
}

// The faces are found without a geometry context, so an alignable volume
// keeps its full portal surfaces. Containers cannot stack alignable volumes,
// so the oversized portal is attached by hand.

BOOST_AUTO_TEST_CASE(AlignableVolumeHasNoFaces) {
  auto bounds = std::make_shared<CylinderVolumeBounds>(rMin, rMax, halfZ);
  auto makeGeometry = [&](std::shared_ptr<TrackingVolume> volume) {
    const auto volumeId = GeometryIdentifier().withVolume(1);
    volume->assignGeometryId(volumeId);
    auto surface = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                        rMax, 3 * halfZ);
    surface->assignGeometryId(volumeId.withBoundary(1));
    volume->addPortal(std::make_shared<Portal>(Direction::OppositeNormal(),
                                               std::move(surface), *volume));
    return std::make_shared<TrackingGeometry>(
        volume, nullptr, GeometryIdentifierHook{}, getDummyLogger(), false);
  };

  auto fixed =
      std::make_shared<TrackingVolume>(Transform3::Identity(), bounds, "fixed");
  auto fixedGeometry = makeGeometry(fixed);
  BOOST_CHECK(fixed->portalFace(*fixed->portals().begin()) != nullptr);

  auto placement = std::make_shared<FixedPlacement>(Transform3::Identity());
  auto alignable =
      std::make_shared<TrackingVolume>(*placement, bounds, "alignable");
  alignable->retainPlacement(placement);
  BOOST_REQUIRE(alignable->isAlignable());
  auto alignableGeometry = makeGeometry(alignable);
  BOOST_CHECK(alignable->portalFace(*alignable->portals().begin()) == nullptr);
}

BOOST_AUTO_TEST_CASE(CuboidStackFaces) {
  auto logger = getDefaultLogger("PortalFace", Logging::INFO);
  Blueprint::Config cfg;
  cfg.envelope = ExtentEnvelope{{
      .x = {20_mm, 20_mm},
      .y = {20_mm, 20_mm},
      .z = {20_mm, 20_mm},
  }};
  Blueprint root{cfg};
  constexpr double half = 50_mm;
  auto& stack = root.addCuboidContainer("Stack", AxisDirection::AxisX);
  auto bounds = std::make_shared<CuboidVolumeBounds>(half, half, half);
  for (int i = 0; i < 2; ++i) {
    stack.addStaticVolume(std::make_unique<TrackingVolume>(
        Transform3{Translation3{Vector3{(2 * i - 1) * half, 0, 0}}}, bounds,
        "box" + std::to_string(i)));
  }
  auto geometry = root.construct({}, gctx, *logger);

  for (int i = 0; i < 2; ++i) {
    const TrackingVolume& volume =
        findVolume(*geometry, "box" + std::to_string(i));
    int merged = 0;
    for (const Portal& portal : volume.portals()) {
      const Surface* face = volume.portalFace(portal);
      // The faces normal to the stacking axis are not merged
      const Vector3 normal = portal.surface().normal(
          gctx, portal.surface().center(gctx), Vector3::UnitX());
      if (std::abs(normal.x()) > 0.5) {
        BOOST_CHECK(face == nullptr);
        continue;
      }
      BOOST_REQUIRE(face != nullptr);
      ++merged;
      const auto& faceBounds =
          dynamic_cast<const RectangleBounds&>(face->bounds());
      BOOST_CHECK_CLOSE(faceBounds.halfLengthX(), half, 1e-9);
      BOOST_CHECK_CLOSE(faceBounds.halfLengthY(), half, 1e-9);
      BOOST_CHECK_SMALL(
          (face->center(gctx) - volume.center(gctx)).norm() - half, 1e-9);
    }
    BOOST_CHECK_EQUAL(merged, 4);
  }
}

// A straight line from inside the first volume reaches the merged outer
// cylinder only beyond the end of that volume

BOOST_AUTO_TEST_CASE(StreamRejectsHitOutsideTheFace) {
  auto logger = getDefaultLogger("PortalFace", Logging::INFO);
  auto geometry = makeCylinderStack(*logger);
  const TrackingVolume& volume = findVolume(*geometry, "child0");
  const Portal& portal = outerCylinder(volume);

  // Reaches r = rMax at z = -2 halfZ + 40 mm, beyond the face at -halfZ
  const NavigationStream::QueryPoint query{
      .position = Vector3{15_mm, 0, -2 * halfZ},
      .direction = Vector3{5_mm, 0, 40_mm}.normalized()};

  NavigationStream merged;
  AppendOnlyNavigationStream{merged}.addPortalCandidate(portal);
  merged.initialize(gctx, query, *logger);
  BOOST_CHECK_EQUAL(merged.candidates().size(), 1u);

  NavigationStream clipped;
  AppendOnlyNavigationStream{clipped, &volume}.addPortalCandidate(portal);
  clipped.initialize(gctx, query, *logger);
  BOOST_CHECK_EQUAL(clipped.candidates().size(), 0u);
}

// The navigator targets the face, but reports the portal surface on arrival

BOOST_AUTO_TEST_CASE(NavigatorReachesPortalSurface) {
  auto logger = getDefaultLogger("PortalFace", Logging::INFO);
  auto geometry = makeCylinderStack(*logger);
  const TrackingVolume& volume = findVolume(*geometry, "child0");
  const Portal& portal = outerCylinder(volume);

  Navigator::Config cfg;
  cfg.trackingGeometry = geometry;
  Navigator navigator(cfg, logger->clone("Navigator"));
  Navigator::Options options(gctx);
  Navigator::State state = navigator.makeState(options);

  Vector3 position{15_mm, 0, -2 * halfZ};
  const Vector3 direction = Vector3::UnitX();
  BOOST_REQUIRE(
      navigator
          .initialize(state, {.position = position, .direction = direction})
          .ok());
  BOOST_REQUIRE_EQUAL(navigator.currentVolume(state), &volume);

  NavigationTarget target = navigator.nextTarget(state, position, direction);
  BOOST_REQUIRE(target.isPortalTarget());
  BOOST_CHECK_EQUAL(&target.portal(), &portal);
  BOOST_CHECK_EQUAL(&target.surface(), volume.portalFace(portal));

  position += target.intersection().pathLength() * direction;
  navigator.handleSurfaceReached(state, position, direction, target.surface());
  BOOST_CHECK_EQUAL(navigator.currentSurface(state), &portal.surface());
  BOOST_CHECK_NE(navigator.currentVolume(state), &volume);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
