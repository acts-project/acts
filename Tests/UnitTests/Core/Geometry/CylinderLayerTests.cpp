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
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/GenericApproachDescriptor.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace Acts;

namespace ActsTests {

namespace {

auto makeApproachTestLayer(const std::shared_ptr<Surface>& sensor,
                           bool cylindrical = true) {
  sensor->assignGeometryId(GeometryIdentifier().withLayer(1).withSensitive(1));
  std::vector<std::shared_ptr<const Surface>> approaches;
  for (double radius : {9., 11.}) {
    auto approach = Surface::makeShared<CylinderSurface>(
        Transform3::Identity(), std::make_shared<CylinderBounds>(radius, 100.));
    approach->assignGeometryId(
        GeometryIdentifier().withLayer(1).withApproach(approaches.size() + 1));
    approach->assignSurfaceMaterial(
        std::make_shared<HomogeneousSurfaceMaterial>());
    approaches.push_back(approach);
  }
  auto array = std::make_unique<SurfaceArray>(sensor);
  auto descriptor =
      std::make_unique<GenericApproachDescriptor>(std::move(approaches));
  std::shared_ptr<Layer> layer;
  if (cylindrical) {
    layer = CylinderLayer::create(Transform3::Identity(),
                                  std::make_shared<CylinderBounds>(10., 100.),
                                  std::move(array), 2., std::move(descriptor));
  } else {
    layer = PlaneLayer::create(Transform3::Identity(),
                               std::make_shared<RectangleBounds>(100., 100.),
                               std::move(array), 2., std::move(descriptor));
  }
  return layer;
}

auto unreachableSensor() {
  return Surface::makeShared<PlaneSurface>(
      Transform3(Translation3(1000., 0., 0.)),
      std::make_shared<RectangleBounds>(1., 1.));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(GeometrySuite)

/// Unit test for creating compliant/non-compliant CylinderLayer object
BOOST_AUTO_TEST_CASE(CylinderLayerConstruction) {
  // default constructor, copy and assignment are all deleted
  // minimally need a Transform3 and a PlanarBounds object (e.g.
  // CylinderBounds) to construct
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
  double radius(0.5), halfz(10.);
  auto pCylinder = std::make_shared<const CylinderBounds>(radius, halfz);
  auto pCylinderLayer =
      CylinderLayer::create(pTransform, pCylinder, nullptr, 1.);
  BOOST_CHECK_EQUAL(pCylinderLayer->layerType(), LayerType::passive);
  // next level: need an array of Surfaces;
  // bounds object, rectangle type
  auto rBounds = std::make_shared<const RectangleBounds>(1., 1.);
  /// Construction
  const std::vector<std::shared_ptr<const Surface>> aSurfaces{
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), rBounds),
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), rBounds)};
  const double thickness(1.0);
  auto pCylinderLayerFromSurfaces =
      CylinderLayer::create(pTransform, pCylinder, nullptr, thickness);
  BOOST_CHECK_EQUAL(pCylinderLayerFromSurfaces->layerType(),
                    LayerType::passive);
  // construct with thickness:
  auto pCylinderLayerWithThickness =
      CylinderLayer::create(pTransform, pCylinder, nullptr, thickness);
  CHECK_CLOSE_REL(pCylinderLayerWithThickness->layerThickness(), thickness,
                  1e-6);
  // with an approach descriptor...
  auto ad(std::make_unique<GenericApproachDescriptor>(aSurfaces));
  auto adPtr = ad.get();
  auto pCylinderLayerWithApproachDescriptor = CylinderLayer::create(
      pTransform, pCylinder, nullptr, thickness, std::move(ad));
  BOOST_CHECK_EQUAL(pCylinderLayerWithApproachDescriptor->approachDescriptor(),
                    adPtr);
  // with the layerType specified...
  auto pCylinderLayerWithLayerType =
      CylinderLayer::create(pTransform, pCylinder, nullptr, thickness,
                            std::move(ad), LayerType::passive);
  BOOST_CHECK_EQUAL(pCylinderLayerWithLayerType->layerType(),
                    LayerType::passive);
}

/// Unit test for testing Layer properties
BOOST_AUTO_TEST_CASE(CylinderLayerProperties) {
  Translation3 translation{0., 1., 2.};
  auto pTransform = Transform3(translation);
  double radius(0.5), halfz(10.);
  auto pCylinder = std::make_shared<const CylinderBounds>(radius, halfz);
  auto pCylinderLayer =
      CylinderLayer::create(pTransform, pCylinder, nullptr, 1.);
  // auto planeSurface = pCylinderLayer->surfaceRepresentation();
  BOOST_CHECK_EQUAL(pCylinderLayer->surfaceRepresentation().name(),
                    std::string("Acts::CylinderSurface"));
}

BOOST_AUTO_TEST_CASE(CylinderLayerSensitiveCylinderKeepsClosestSelection) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto sensor = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), std::make_shared<CylinderBounds>(10., 90.));
  auto layer = makeApproachTestLayer(sensor);
  // Resolve the material approaches too, so side selection is active.
  NavigationOptions<Surface> options;
  const auto targets =
      layer->compatibleSurfaces(gctx, Vector3(9.5, 0., 97.5),
                                Vector3(-0.2, 0., -1.).normalized(), options);
  const auto found = std::ranges::find_if(targets, [&](const auto& target) {
    return &target.surface() == sensor.get();
  });
  // The approach has root0 ahead, but only the sensor's far root is in bounds.
  BOOST_REQUIRE(found != targets.end());
  BOOST_CHECK_EQUAL(found->intersectionIndex(), 1u);
  CHECK_CLOSE_ABS(found->position().x(), -10., 1e-9);
}

BOOST_DATA_TEST_CASE(CylinderLayerApproachSideIncludesStartSurface,
                     boost::unit_test::data::make({false, true}), onApproach) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto layer = makeApproachTestLayer(unreachableSensor());
  const Surface* inner =
      layer->approachDescriptor()->containedSurfaces().front();
  NavigationOptions<Surface> options;
  options.resolvePassive = true;
  options.nearLimit = s_onSurfaceTolerance;
  options.startObject = onApproach ? inner : nullptr;
  // The representing cylinder's near root is out of bounds as well, so
  // passive collection must not substitute its far root either.
  const Vector3 position =
      onApproach ? Vector3(9., 0., 95.01) : Vector3(9.5, 0., 97.51);
  const auto targets = layer->compatibleSurfaces(
      gctx, position, Vector3(-0.2, 0., -1.).normalized(), options);
  BOOST_REQUIRE_EQUAL(targets.size(), onApproach ? 0u : 1u);
  if (!onApproach) {
    BOOST_CHECK_EQUAL(&targets.front().surface(), inner);
    BOOST_CHECK_EQUAL(targets.front().intersectionIndex(), 0u);
  }
}

BOOST_AUTO_TEST_CASE(CylinderLayerApproachSideUsesExtendedBounds) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto layer = makeApproachTestLayer(unreachableSensor());
  const Surface* inner =
      layer->approachDescriptor()->containedSurfaces().front();
  NavigationOptions<Surface> options;
  options.resolveSensitive = false;
  options.extendedSurfaces.emplace_back(inner->geometryId(),
                                        BoundaryTolerance::Infinite());
  const auto targets =
      layer->compatibleSurfaces(gctx, Vector3(9.5, 0., 110.),
                                Vector3(-0.2, 0., -1.).normalized(), options);
  BOOST_REQUIRE_EQUAL(targets.size(), 1u);
  BOOST_CHECK_EQUAL(&targets.front().surface(), inner);
  BOOST_CHECK_EQUAL(targets.front().intersectionIndex(), 0u);
  BOOST_CHECK(targets.front().boundaryTolerance().isInfinite());
}

BOOST_AUTO_TEST_CASE(
    CylinderLayerWithoutFirstApproachRootKeepsClosestSelection) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto layer = makeApproachTestLayer(unreachableSensor());
  NavigationOptions<Surface> options;
  options.resolveSensitive = false;
  const auto outward = layer->compatibleSurfaces(gctx, Vector3(9.5, 0., 0.),
                                                 Vector3::UnitX(), options);
  BOOST_REQUIRE_EQUAL(outward.size(), 1u);
  BOOST_CHECK_EQUAL(outward.front().intersectionIndex(), 1u);
  CHECK_CLOSE_ABS(outward.front().position().x(), 11., 1e-9);
  // Inner cylinder is missed; the outer cylinder's closest root is behind.
  // Having no eligible root0 must not force selection of its far root1.
  BOOST_CHECK(layer
                  ->compatibleSurfaces(gctx, Vector3(-5., 9.5, 0.),
                                       Vector3::UnitX(), options)
                  .empty());
}

BOOST_AUTO_TEST_CASE(PlaneLayerCylinderApproachesKeepClosestSelection) {
  const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
  auto layer = makeApproachTestLayer(unreachableSensor(), false);
  NavigationOptions<Surface> options;
  options.resolveSensitive = false;
  const auto targets =
      layer->compatibleSurfaces(gctx, Vector3(9.5, 0., 97.5),
                                Vector3(-0.2, 0., -1.).normalized(), options);
  BOOST_REQUIRE_EQUAL(targets.size(), 2u);
  BOOST_CHECK_EQUAL(targets[0].intersectionIndex(), 0u);
  BOOST_CHECK_EQUAL(targets[1].intersectionIndex(), 1u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
