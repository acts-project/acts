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
#include "Acts/Geometry/GeometryIdentifierBlueprintNode.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/LinearVolumeLookup.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/PadBlueprintNode.hpp"
#include "Acts/Geometry/PortalDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/VolumePlacementBase.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsTests/CommonHelpers/DetectorElementStub.hpp"

#include <map>
#include <memory>
#include <optional>
#include <string>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace {
const auto gctx = GeometryContext::dangerouslyDefaultConstruct();

// Deliberately keeps a derived surface private: a post-construction visitor
// cannot repair its placement. It must be created in the final frame.
class CachedSurfacePolicy final : public INavigationPolicy {
 public:
  CachedSurfacePolicy(const GeometryContext& context,
                      const TrackingVolume& volume, const Logger& /*logger*/)
      : m_surface(Surface::makeShared<PlaneSurface>(
            volume.localToGlobalTransform(context),
            std::make_shared<RectangleBounds>(5_mm, 5_mm))) {}

  void connect(NavigationDelegate& delegate) const override {
    connectDefault<CachedSurfacePolicy>(delegate);
  }

  void initializeCandidates(const GeometryContext&, const NavigationArguments&,
                            NavigationPolicyState&,
                            AppendOnlyNavigationStream& stream,
                            const Logger&) const {
    stream.addSurfaceCandidate(*m_surface, BoundaryTolerance::None());
  }

 private:
  std::shared_ptr<PlaneSurface> m_surface;
};

std::unique_ptr<TrackingGeometry> telescope(
    std::optional<Transform3> placement) {
  const ExtentEnvelope envelope{
      {.x = {1_mm, 1_mm}, .y = {1_mm, 1_mm}, .z = {1_mm, 1_mm}}};
  Blueprint root{{.envelope = envelope, .globalTransform = placement}};
  auto& ids = root.addMaterial("Material")
                  .addPortalDesignator("Portals")
                  .withGeometryIdentifier();
  ids.setDirectChildVolumeIdTo(3);
  auto pad = std::make_shared<PadBlueprintNode>("Envelope", envelope);
  ids.addChild(pad);
  const auto lookupFactory = [placement](const GeometryContext& context,
                                         const TrackingVolume& volume,
                                         const Logger& logger) {
    BOOST_CHECK(volume.localToGlobalTransform(context).linear().isApprox(
        placement.value_or(Transform3::Identity()).linear()));
    return std::make_unique<LinearVolumeLookup>(
        volume, LinearVolumeLookup::Config{}, logger);
  };
  root.setVolumeLookupFactory(lookupFactory);
  pad->setVolumeLookupFactory(lookupFactory);
  auto& outer = pad->addCuboidContainer("Outer", AxisDirection::AxisZ);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
  auto& inner = outer.addCuboidContainer("Inner", AxisDirection::AxisZ);
  inner.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
  for (int i = 0; i < 3; ++i) {
    BlueprintNode& parent = i < 2 ? static_cast<BlueprintNode&>(inner) : outer;
    auto surface = Surface::makeShared<PlaneSurface>(
        Transform3{Translation3{0., 0., i * 20_mm}},
        std::make_shared<RectangleBounds>(10_mm, 20_mm));
    surface->assignIsSensitive(true);
    parent.addLayer("Plane" + std::to_string(i))
        .setLayerType(LayerBlueprintNode::LayerType::Plane)
        .setSurfaces({surface})
        .setEnvelope(envelope);
  }
  BlueprintOptions options;
  options.defaultNavigationPolicyFactory =
      NavigationPolicyFactory{}.add<CachedSurfacePolicy>().asUniquePtr();
  return root.construct(options, gctx);
}

std::map<std::string, const TrackingVolume*> volumes(
    const TrackingGeometry& geometry) {
  std::map<std::string, const TrackingVolume*> result;
  geometry.visitVolumes([&](const TrackingVolume* volume) {
    result.emplace(volume->volumeName(), volume);
  });
  return result;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(BlueprintPlacement)

BOOST_AUTO_TEST_CASE(RotatedTelescope) {
  const Transform3 placement =
      Translation3{100_mm, -70_mm, 40_mm} *
      AngleAxis3{120_degree, Vector3{1., 2., 3.}.normalized()};
  auto nominal = telescope(std::nullopt);
  auto placed = telescope(placement);
  auto referenceVolumes = volumes(*nominal);
  auto placedVolumes = volumes(*placed);
  BOOST_REQUIRE_EQUAL(referenceVolumes.size(), placedVolumes.size());
  BOOST_CHECK(referenceVolumes.contains("Inner::Gap1"));
  BOOST_CHECK(referenceVolumes.contains("Outer::Gap1"));
  BOOST_REQUIRE(placedVolumes.at("World")->volumeLookup() != nullptr);
  BOOST_REQUIRE(placedVolumes.at("Envelope")->volumeLookup() != nullptr);
  for (double z : {0_mm, 10_mm, 20_mm, 30_mm, 40_mm}) {
    const auto* expected =
        nominal->lowestTrackingVolume(gctx, Vector3{0., 0., z});
    const auto* actual =
        placed->lowestTrackingVolume(gctx, placement * Vector3{0., 0., z});
    BOOST_REQUIRE(expected != nullptr);
    BOOST_REQUIRE(actual != nullptr);
    BOOST_CHECK_EQUAL(actual->volumeName(), expected->volumeName());
  }
  for (const auto& [name, reference] : referenceVolumes) {
    const auto& actual = *placedVolumes.at(name);
    BOOST_CHECK(actual.localToGlobalTransform(gctx).isApprox(
        placement * reference->localToGlobalTransform(gctx)));
    BOOST_CHECK_EQUAL(actual.geometryId(), reference->geometryId());
    BOOST_REQUIRE_EQUAL(actual.surfaces().size(), reference->surfaces().size());
    for (std::size_t i = 0; i < actual.surfaces().size(); ++i) {
      BOOST_CHECK(actual.surfaces()[i].localToGlobalTransform(gctx).isApprox(
          placement * reference->surfaces()[i].localToGlobalTransform(gctx)));
      BOOST_CHECK(actual.surfaces()[i].isSensitive());
    }
    BOOST_REQUIRE_EQUAL(actual.portals().size(), reference->portals().size());
    for (std::size_t i = 0; i < actual.portals().size(); ++i) {
      const auto& before = reference->portals()[i];
      const auto& after = actual.portals()[i];
      BOOST_CHECK(after.surface().localToGlobalTransform(gctx).isApprox(
          placement * before.surface().localToGlobalTransform(gctx)));
      const Vector3 point = before.surface().localToGlobalTransform(gctx) *
                            Vector3{0.123_mm, 0.234_mm, 0.};
      const Vector3 normal =
          before.surface().normal(gctx, point, Vector3::UnitZ());
      for (double sign : {-1., 1.}) {
        auto expected =
            before.resolveVolume(gctx, point, sign * normal).value();
        auto resolved = after
                            .resolveVolume(gctx, placement * point,
                                           placement.linear() * (sign * normal))
                            .value();
        BOOST_CHECK_EQUAL(resolved == nullptr, expected == nullptr);
        if (expected != nullptr && resolved != nullptr) {
          BOOST_CHECK_EQUAL(resolved->volumeName(), expected->volumeName());
        }
      }
    }
    // Navigate against the policy's private surface using a transformed ray.
    NavigationArguments args{
        .position =
            actual.localToGlobalTransform(gctx) * Vector3{0., 0., -0.5_mm},
        .direction =
            actual.localToGlobalTransform(gctx).linear() * Vector3::UnitZ()};
    NavigationPolicyStateManager manager;
    const auto* policy = actual.navigationPolicy();
    BOOST_REQUIRE(policy != nullptr);
    policy->createState(gctx, args, manager, getDummyLogger());
    auto state = manager.currentState();
    NavigationDelegate delegate;
    policy->connect(delegate);
    NavigationStream stream;
    AppendOnlyNavigationStream append{stream};
    delegate(gctx, args, state, append, getDummyLogger());
    BOOST_REQUIRE(stream.initialize(
        gctx, {.position = args.position, .direction = args.direction},
        getDummyLogger()));
    BOOST_CHECK_SMALL(
        stream.currentCandidate().intersection().pathLength() - 0.5_mm, 1e-9);
  }
}

BOOST_AUTO_TEST_CASE(SharedSurfacePlacedOnce) {
  const Transform3 placement{Translation3{10_mm, 20_mm, 30_mm}};
  auto surface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(1_mm, 1_mm));
  auto parent = std::make_unique<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(10_mm, 10_mm, 10_mm), "Parent");
  auto child = std::make_unique<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(2_mm, 2_mm, 2_mm), "Child");
  parent->addSurface(surface);
  child->addSurface(surface);
  parent->addVolume(std::move(child));
  Blueprint root{{.globalTransform = placement}};
  root.addStaticVolume(std::move(parent));
  auto geometry = root.construct({}, gctx);
  BOOST_CHECK(surface->localToGlobalTransform(gctx).isApprox(placement));
  BOOST_CHECK(
      volumes(*geometry).at("Child")->localToGlobalTransform(gctx).isApprox(
          placement));
}

BOOST_AUTO_TEST_CASE(RejectExternalPlacementBeforeMutation) {
  ActsTests::DetectorElementStub element{
      Transform3::Identity(), std::make_shared<RectangleBounds>(1_mm, 1_mm),
      0.1_mm};
  auto fixed = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(1_mm, 1_mm));
  auto volume = std::make_unique<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(10_mm, 10_mm, 10_mm), "Volume");
  auto* volumePtr = volume.get();
  volume->addSurface(fixed);
  volume->addSurface(element.surface().getSharedPtr());
  const Transform3 placement{Translation3{10_mm, 20_mm, 30_mm}};
  Blueprint root{{.globalTransform = placement}};
  root.addStaticVolume(std::move(volume));
  BOOST_CHECK_THROW(root.construct({}, gctx), std::invalid_argument);
  BOOST_CHECK(
      volumePtr->localToGlobalTransform(gctx).isApprox(Transform3::Identity()));
  BOOST_CHECK(
      fixed->localToGlobalTransform(gctx).isApprox(Transform3::Identity()));
  BOOST_CHECK_THROW(element.surface().setTransform(placement),
                    std::runtime_error);
  BOOST_CHECK(element.surface().isSensitive());
}

BOOST_AUTO_TEST_CASE(RejectExternalVolumePlacement) {
  class Placement final : public VolumePlacementBase {
   public:
    const Transform3& localToGlobalTransform(
        const GeometryContext&) const override {
      return transform;
    }
    const Transform3& globalToLocalTransform(
        const GeometryContext&) const override {
      return transform;
    }
    const Transform3& portalLocalToGlobal(const GeometryContext&,
                                          std::size_t) const override {
      throw std::logic_error(
          "Portals must not be created before placement validation");
    }
    Transform3 transform = Transform3::Identity();
  } external;
  Blueprint root{
      {.globalTransform = Transform3{Translation3{10_mm, 20_mm, 30_mm}}}};
  root.addStaticVolume(std::make_unique<TrackingVolume>(
      external, std::make_shared<CuboidVolumeBounds>(10_mm, 10_mm, 10_mm),
      "Aligned"));
  BOOST_CHECK_THROW(root.construct({}, gctx), std::invalid_argument);
  BOOST_CHECK(external.transform.isApprox(Transform3::Identity()));
  BOOST_CHECK_EQUAL(external.nPortalPlacements(), 0u);
}

BOOST_AUTO_TEST_CASE(RejectPreexistingPolicy) {
  auto volume = std::make_unique<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(10_mm, 10_mm, 10_mm), "Volume");
  auto* original = volume.get();
  volume->setNavigationPolicy(
      std::make_unique<CachedSurfacePolicy>(gctx, *volume, getDummyLogger()));
  Blueprint root{
      {.globalTransform = Transform3{Translation3{10_mm, 20_mm, 30_mm}}}};
  root.addStaticVolume(std::move(volume));
  BOOST_CHECK_THROW(root.construct({}, gctx), std::invalid_argument);
  BOOST_CHECK(
      original->localToGlobalTransform(gctx).isApprox(Transform3::Identity()));
}

BOOST_AUTO_TEST_CASE(RejectPreexistingVolumeLookup) {
  auto volume = std::make_unique<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(10_mm, 10_mm, 10_mm), "Volume");
  auto* original = volume.get();
  volume->setVolumeLookup(std::make_unique<LinearVolumeLookup>(
      *volume, LinearVolumeLookup::Config{}));
  Blueprint root{
      {.globalTransform = Transform3{Translation3{10_mm, 20_mm, 30_mm}}}};
  root.addStaticVolume(std::move(volume));
  BOOST_CHECK_THROW(root.construct({}, gctx), std::invalid_argument);
  BOOST_CHECK(
      original->localToGlobalTransform(gctx).isApprox(Transform3::Identity()));
}

BOOST_AUTO_TEST_CASE(RotatedCylinderStack) {
  const Transform3 placement = Translation3{30_mm, 50_mm, -20_mm} *
                               AngleAxis3{90_degree, Vector3::UnitY()};
  auto make = [&](std::optional<Transform3> transform) {
    Blueprint root{
        {.envelope = ExtentEnvelope{{.z = {1_mm, 1_mm}, .r = {1_mm, 1_mm}}},
         .globalTransform = transform}};
    auto& stack = root.addCylinderContainer("Cylinders", AxisDirection::AxisZ);
    stack.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
    for (int i = 0; i < 2; ++i) {
      stack.addStaticVolume(
          Transform3{Translation3{0., 0., i * 30_mm}},
          std::make_shared<CylinderVolumeBounds>(5_mm, 10_mm, 10_mm),
          "Cylinder" + std::to_string(i));
    }
    return root.construct({}, gctx);
  };
  auto nominal = make(std::nullopt);
  auto placed = make(placement);
  auto actual = volumes(*placed);
  BOOST_CHECK(actual.contains("Cylinders::Gap1"));
  for (const auto& [name, reference] : volumes(*nominal)) {
    auto* rotated = actual.at(name);
    BOOST_CHECK(rotated->localToGlobalTransform(gctx).isApprox(
        placement * reference->localToGlobalTransform(gctx)));
    BOOST_REQUIRE_EQUAL(rotated->portals().size(), reference->portals().size());
    for (std::size_t i = 0; i < rotated->portals().size(); ++i) {
      BOOST_CHECK(
          rotated->portals()[i].surface().localToGlobalTransform(gctx).isApprox(
              placement *
              reference->portals()[i].surface().localToGlobalTransform(gctx)));
    }
  }
}

BOOST_AUTO_TEST_CASE(RejectNonRigidPlacement) {
  Transform3 scale = Transform3::Identity();
  scale.linear() *= 2.;
  Blueprint root{{.globalTransform = scale}};
  root.addStaticVolume(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(10_mm, 10_mm, 10_mm), "Volume");
  BOOST_CHECK_THROW(root.construct({}, gctx), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
