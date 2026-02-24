// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Geometry/CompositePortalLink.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"
#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"

#include <functional>
#include <fstream>
#include <memory>
#include <utility>
#include <vector>

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(TrackingGeometryJsonConverterRoundTrip) {
  using namespace Acts;

  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

  auto root = std::make_shared<TrackingVolume>(
      Transform3::Identity(), std::make_shared<CuboidVolumeBounds>(5., 5., 5.),
      "root");
  root->assignGeometryId(GeometryIdentifier{}.withVolume(1u));
  TrackingVolume* rootPtr = root.get();

  Transform3 childTransform = Transform3::Identity();
  childTransform.pretranslate(Vector3{1., 0., 0.});
  auto child = std::make_unique<TrackingVolume>(
      childTransform, std::make_shared<CylinderVolumeBounds>(0.5, 1.0, 2.0),
      "child");
  child->assignGeometryId(GeometryIdentifier{}.withVolume(2u));
  TrackingVolume* childPtr = child.get();
  root->addVolume(std::move(child));

  auto trivialBounds = std::make_shared<const RectangleBounds>(1., 1.);
  auto trivialSurface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), trivialBounds);
  trivialSurface->assignGeometryId(
      GeometryIdentifier{}.withVolume(1u).withExtra(11u));
  auto trivialLink =
      std::make_unique<TrivialPortalLink>(trivialSurface, *childPtr);
  root->addPortal(
      std::make_shared<Portal>(Direction::AlongNormal(), std::move(trivialLink)));

  auto compositeBounds = std::make_shared<const RectangleBounds>(1., 1.);
  Transform3 compositeTransformA = Transform3::Identity();
  compositeTransformA.pretranslate(Vector3{-1., 0., 0.});
  Transform3 compositeTransformB = Transform3::Identity();
  compositeTransformB.pretranslate(Vector3{1., 0., 0.});
  auto compositeSurfaceA = Surface::makeShared<PlaneSurface>(
      compositeTransformA, compositeBounds);
  auto compositeSurfaceB = Surface::makeShared<PlaneSurface>(
      compositeTransformB, compositeBounds);
  compositeSurfaceA->assignGeometryId(
      GeometryIdentifier{}.withVolume(1u).withExtra(12u));
  compositeSurfaceB->assignGeometryId(
      GeometryIdentifier{}.withVolume(1u).withExtra(13u));
  auto compositeLinkA =
      std::make_unique<TrivialPortalLink>(compositeSurfaceA, *childPtr);
  auto compositeLinkB =
      std::make_unique<TrivialPortalLink>(compositeSurfaceB, *rootPtr);
  auto compositeLink = std::make_unique<CompositePortalLink>(
      std::move(compositeLinkA), std::move(compositeLinkB),
      AxisDirection::AxisX);
  root->addPortal(std::make_shared<Portal>(Direction::AlongNormal(),
                                           std::move(compositeLink)));

  auto gridBounds = std::make_shared<const RectangleBounds>(2., 1.);
  auto gridSurface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), gridBounds);
  gridSurface->assignGeometryId(
      GeometryIdentifier{}.withVolume(1u).withExtra(14u));
  auto gridLink = GridPortalLink::make(
      gridSurface, AxisDirection::AxisX, Axis{AxisBound, -2., 2., 2});
  AnyGridView<const TrackingVolume*> gridView(gridLink->grid());
  gridView.atLocalBins({0u}) = rootPtr;
  gridView.atLocalBins({1u}) = childPtr;
  gridView.atLocalBins({2u}) = rootPtr;
  gridView.atLocalBins({3u}) = childPtr;

  std::vector<TrivialPortalLink> artifacts;
  auto artifactSurface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), gridBounds);
  artifactSurface->assignGeometryId(
      GeometryIdentifier{}.withVolume(1u).withExtra(15u));
  artifacts.emplace_back(artifactSurface, *childPtr);
  gridLink->setArtifactPortalLinks(std::move(artifacts));

  root->addPortal(
      std::make_shared<Portal>(Direction::AlongNormal(), std::move(gridLink)));

  auto sharedPortalBounds = std::make_shared<const RectangleBounds>(0.75, 0.75);
  auto sharedPortalSurface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), sharedPortalBounds);
  sharedPortalSurface->assignGeometryId(
      GeometryIdentifier{}.withVolume(1u).withExtra(16u));

  auto sharedPortal = std::make_shared<Portal>(
      gctx, std::make_unique<TrivialPortalLink>(sharedPortalSurface, *childPtr),
      std::make_unique<TrivialPortalLink>(sharedPortalSurface, *rootPtr));
  root->addPortal(sharedPortal);
  childPtr->addPortal(sharedPortal);

  TrackingGeometryJsonConverter converter;
  nlohmann::json encoded = converter.toJson(gctx, *root);
  TemporaryDirectory tmpDir{};
  auto jsonPath = tmpDir.path() / "tracking_geometry_roundtrip.json";

  {
    std::ofstream out(jsonPath);
    BOOST_REQUIRE(out.good());
    out << encoded.dump(2);
  }

  nlohmann::json encodedFromFile;
  {
    std::ifstream in(jsonPath);
    BOOST_REQUIRE(in.good());
    in >> encodedFromFile;
  }

  auto decodedRoot = converter.trackingVolumeFromJson(gctx, encodedFromFile);

  BOOST_REQUIRE(decodedRoot != nullptr);
  BOOST_CHECK_EQUAL(decodedRoot->volumeName(), "root");
  BOOST_CHECK_EQUAL(decodedRoot->volumeBounds().type(), VolumeBounds::eCuboid);

  std::vector<TrackingVolume*> decodedChildren;
  for (auto& decodedChild : decodedRoot->volumes()) {
    decodedChildren.push_back(&decodedChild);
  }
  BOOST_REQUIRE_EQUAL(decodedChildren.size(), 1u);
  BOOST_CHECK_EQUAL(decodedChildren.front()->volumeName(), "child");
  BOOST_CHECK_EQUAL(decodedChildren.front()->volumeBounds().type(),
                    VolumeBounds::eCylinder);

  std::vector<Portal*> decodedPortals;
  for (auto& portal : decodedRoot->portals()) {
    decodedPortals.push_back(&portal);
  }
  BOOST_REQUIRE_EQUAL(decodedPortals.size(), 4u);

  const auto* decodedTrivial = dynamic_cast<const TrivialPortalLink*>(
      decodedPortals.at(0)->getLink(Direction::AlongNormal()));
  BOOST_REQUIRE(decodedTrivial != nullptr);
  BOOST_CHECK_EQUAL(decodedTrivial->volume().volumeName(), "child");

  const auto* decodedComposite = dynamic_cast<const CompositePortalLink*>(
      decodedPortals.at(1)->getLink(Direction::AlongNormal()));
  BOOST_REQUIRE(decodedComposite != nullptr);
  BOOST_CHECK_EQUAL(decodedComposite->size(), 2u);
  BOOST_CHECK_EQUAL(decodedComposite->direction(), AxisDirection::AxisX);

  const auto* decodedGrid = dynamic_cast<const GridPortalLink*>(
      decodedPortals.at(2)->getLink(Direction::AlongNormal()));
  BOOST_REQUIRE(decodedGrid != nullptr);
  BOOST_CHECK_EQUAL(decodedGrid->dim(), 1u);
  BOOST_CHECK_EQUAL(decodedGrid->artifactPortalLinks().size(), 1u);

  AnyGridConstView<const TrackingVolume*> decodedGridView(decodedGrid->grid());
  BOOST_REQUIRE(decodedGridView.atLocalBins({0u}) != nullptr);
  BOOST_REQUIRE(decodedGridView.atLocalBins({1u}) != nullptr);
  BOOST_REQUIRE(decodedGridView.atLocalBins({2u}) != nullptr);
  BOOST_REQUIRE(decodedGridView.atLocalBins({3u}) != nullptr);

  BOOST_CHECK_EQUAL(decodedGridView.atLocalBins({0u})->volumeName(), "root");
  BOOST_CHECK_EQUAL(decodedGridView.atLocalBins({1u})->volumeName(), "child");
  BOOST_CHECK_EQUAL(decodedGridView.atLocalBins({2u})->volumeName(), "root");
  BOOST_CHECK_EQUAL(decodedGridView.atLocalBins({3u})->volumeName(), "child");

  std::vector<Portal*> decodedChildPortals;
  for (auto& portal : decodedChildren.front()->portals()) {
    decodedChildPortals.push_back(&portal);
  }
  BOOST_REQUIRE_EQUAL(decodedChildPortals.size(), 1u);

  bool sharedPortalPreserved = false;
  for (Portal* rootPortal : decodedPortals) {
    if (rootPortal == decodedChildPortals.front()) {
      sharedPortalPreserved = true;
      break;
    }
  }
  BOOST_CHECK(sharedPortalPreserved);

  auto decodedGeometry =
      converter.trackingGeometryFromJson(gctx, encodedFromFile);
  BOOST_REQUIRE(decodedGeometry != nullptr);
  BOOST_REQUIRE(decodedGeometry->highestTrackingVolume() != nullptr);
  BOOST_CHECK_EQUAL(decodedGeometry->highestTrackingVolume()->volumeName(),
                    "root");
}

BOOST_AUTO_TEST_CASE(TrackingGeometryJsonConverterRoundTripGen3Cylindrical) {
  using namespace Acts;

  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

  CylindricalTrackingGeometry cylindricalGeometryBuilder(gctx, true);
  auto sourceGeometry = cylindricalGeometryBuilder();

  BOOST_REQUIRE(sourceGeometry != nullptr);
  BOOST_REQUIRE(sourceGeometry->highestTrackingVolume() != nullptr);
  BOOST_CHECK(sourceGeometry->geometryVersion() ==
              TrackingGeometry::GeometryVersion::Gen3);

  auto countVolumesAndPortals = [](const TrackingVolume& world) {
    std::size_t volumeCount = 0u;
    std::size_t portalCount = 0u;

    std::function<void(const TrackingVolume&)> traverse =
        [&](const TrackingVolume& volume) {
          ++volumeCount;
          for (const auto& portal : volume.portals()) {
            static_cast<void>(portal);
            ++portalCount;
          }
          for (const auto& child : volume.volumes()) {
            traverse(child);
          }
        };

    traverse(world);
    return std::pair{volumeCount, portalCount};
  };

  const auto [sourceVolumeCount, sourcePortalCount] =
      countVolumesAndPortals(*sourceGeometry->highestTrackingVolume());
  BOOST_CHECK_GT(sourceVolumeCount, 0u);
  BOOST_CHECK_GT(sourcePortalCount, 0u);

  TrackingGeometryJsonConverter converter;
  nlohmann::json encoded = converter.toJson(gctx, *sourceGeometry);

  TemporaryDirectory tmpDir{};
  auto jsonPath = tmpDir.path() / "tracking_geometry_gen3_roundtrip.json";

  {
    std::ofstream out(jsonPath);
    BOOST_REQUIRE(out.good());
    out << encoded.dump(2);
  }

  nlohmann::json encodedFromFile;
  {
    std::ifstream in(jsonPath);
    BOOST_REQUIRE(in.good());
    in >> encodedFromFile;
  }

  auto decodedGeometry =
      converter.trackingGeometryFromJson(gctx, encodedFromFile);
  BOOST_REQUIRE(decodedGeometry != nullptr);
  BOOST_REQUIRE(decodedGeometry->highestTrackingVolume() != nullptr);

  const auto [decodedVolumeCount, decodedPortalCount] =
      countVolumesAndPortals(*decodedGeometry->highestTrackingVolume());

  BOOST_CHECK_EQUAL(decodedVolumeCount, sourceVolumeCount);
  BOOST_CHECK_EQUAL(decodedPortalCount, sourcePortalCount);
  BOOST_CHECK_EQUAL(decodedGeometry->highestTrackingVolume()->volumeName(),
                    sourceGeometry->highestTrackingVolume()->volumeName());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
