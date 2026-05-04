// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/tree/test_unit.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/CompositePortalLink.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/LayerBlueprintNode.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"

#include <cstddef>
#include <fstream>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace ActsTests {

using namespace Acts;
using namespace Experimental;
using namespace UnitLiterals;
using enum CylinderVolumeBounds::Face;
using enum AxisDirection;

using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

template <typename navigator_t = Navigator,
          typename stepper_t = ConstantFieldStepper>
struct StateCollector {
  struct this_result {
    std::vector<typename navigator_t::State> navigation;
    std::vector<typename stepper_t::State> stepping;
  };

  using result_type = this_result;

  template <typename propagator_state_t>
  Result<void> act(propagator_state_t& state, const stepper_t& /*stepper*/,
                   const navigator_t& /*navigator*/, result_type& result,
                   const Logger& /*logger*/) const {
    result.navigation.push_back(state.navigation);
    result.stepping.push_back(state.stepping);
    return {};
  }
};

void checkHierarchy(const GeometryContext& gctx,
                    const TrackingVolume::VolumeRange& volsA,
                    const TrackingVolume::VolumeRange& volsB) {
  auto checkSurfaces = [&](const auto& surfA, const auto& surfB) {
    BOOST_CHECK_EQUAL(surfA.type(), surfB.type());
    BOOST_CHECK_EQUAL(surfA.geometryId(), surfB.geometryId());
    BOOST_CHECK_EQUAL(surfA.bounds(), surfB.bounds());

    BOOST_CHECK_LT((surfA.localToGlobalTransform(gctx).matrix() -
                    surfB.localToGlobalTransform(gctx).matrix())
                       .norm(),
                   1e-10);
    const auto* matA = surfA.surfaceMaterial();
    const auto* matB = surfB.surfaceMaterial();
    if (matA != nullptr) {
      BOOST_REQUIRE_NE(matB, nullptr);
      BOOST_CHECK_EQUAL(surfA.surfaceMaterial()->materialSlab(Vector2{0, 0}),
                        surfB.surfaceMaterial()->materialSlab(Vector2{0, 0}));
    }
  };
  std::function<void(const Acts::PortalLinkBase* linkA,
                     const Acts::PortalLinkBase* linkB)>
      checkLinks = [&](const auto* linkA, const auto* linkB) {
        const auto* trivialA =
            dynamic_cast<const Acts::TrivialPortalLink*>(linkA);
        const auto* trivialB =
            dynamic_cast<const Acts::TrivialPortalLink*>(linkB);
        if (trivialA != nullptr) {
          BOOST_REQUIRE_NE(trivialB, nullptr);
          BOOST_CHECK(trivialA->volume() == trivialB->volume());
          checkSurfaces(trivialA->surface(), trivialB->surface());
        }

        const auto* gridA = dynamic_cast<const Acts::GridPortalLink*>(linkA);
        const auto* gridB = dynamic_cast<const Acts::GridPortalLink*>(linkB);
        if (gridA != nullptr) {
          BOOST_REQUIRE_NE(gridB, nullptr);
          BOOST_CHECK_EQUAL(gridA->dim(), gridB->dim());
          BOOST_CHECK_EQUAL(gridA->direction(), gridB->direction());
          BOOST_CHECK(gridA->grid() == gridB->grid());
          checkSurfaces(gridA->surface(), gridB->surface());

          const auto& childrenA = gridA->artifactPortalLinks();
          const auto& childrenB = gridA->artifactPortalLinks();
          for (std::size_t i = 0; i < childrenA.size(); i++) {
            const auto& childA = *(childrenA.begin() + i);
            const auto& childB = *(childrenB.begin() + i);
            checkLinks(&childA, &childB);
          }
        }

        const auto* compositeA =
            dynamic_cast<const Acts::CompositePortalLink*>(linkA);
        const auto* compositeB =
            dynamic_cast<const Acts::CompositePortalLink*>(linkB);
        if (compositeA != nullptr) {
          BOOST_REQUIRE_NE(compositeB, nullptr);
          const auto& childrenA = compositeA->links();
          const auto& childrenB = compositeB->links();
          for (std::size_t i = 0; i < childrenA.size(); i++) {
            const auto& childA = childrenA.at(i);
            const auto& childB = childrenB.at(i);
            checkLinks(&childA, &childB);
          }
        }
      };

  BOOST_CHECK_EQUAL(volsA.size(), volsB.size());
  for (std::size_t i = 0; i < volsA.size(); i++) {
    const auto& volA = volsA.at(i);
    const auto& volB = volsB.at(i);
    BOOST_CHECK(volA == volB);

    const auto& surfsA = volA.surfaces();
    const auto& surfsB = volB.surfaces();
    BOOST_CHECK_EQUAL(surfsA.size(), surfsB.size());
    for (std::size_t j = 0; j < surfsA.size(); j++) {
      checkSurfaces(surfsA.at(i), surfsB.at(i));
    }

    const auto& portsA = volA.portals();
    const auto& portsB = volB.portals();
    BOOST_CHECK_EQUAL(portsA.size(), portsB.size());
    for (std::size_t j = 0; j < portsA.size(); j++) {
      const auto& portA = portsA.at(j);
      const auto& portB = portsB.at(j);

      const auto& surfA = portA.surface();
      const auto& surfB = portB.surface();
      checkSurfaces(surfA, surfB);

      const auto* alongA = portA.getLink(Acts::Direction::AlongNormal());
      const auto* alongB = portB.getLink(Acts::Direction::AlongNormal());
      checkLinks(alongA, alongB);

      const auto* oppositeA = portA.getLink(Acts::Direction::OppositeNormal());
      const auto* oppositeB = portB.getLink(Acts::Direction::OppositeNormal());
      checkLinks(oppositeA, oppositeB);
    }
    checkHierarchy(gctx, volA.volumes(), volB.volumes());
  }
}

auto logger = getDefaultLogger("UnitTests", Logging::INFO);

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(TrackingGeometryJsonConverterRoundTrip) {
  using namespace Acts;

  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

  TryAllNavigationPolicy::Config tryAllConfig;
  tryAllConfig.portals = true;
  tryAllConfig.sensitives = false;
  std::unique_ptr<NavigationPolicyFactory> navPolicyFactory =
      NavigationPolicyFactory{}
          .add<TryAllNavigationPolicy>(tryAllConfig)
          .asUniquePtr();

  auto root = std::make_shared<TrackingVolume>(
      Transform3::Identity(), std::make_shared<CuboidVolumeBounds>(5., 5., 5.),
      "root");
  root->assignGeometryId(GeometryIdentifier{}.withVolume(1u));
  TrackingVolume* rootPtr = root.get();
  rootPtr->setNavigationPolicy(
      navPolicyFactory->build(gctx, *rootPtr, *logger));

  Transform3 childTransform = Transform3::Identity();
  childTransform.pretranslate(Vector3{1., 0., 0.});
  auto child = std::make_unique<TrackingVolume>(
      childTransform, std::make_shared<CylinderVolumeBounds>(0.5, 1.0, 2.0),
      "child");
  child->assignGeometryId(GeometryIdentifier{}.withVolume(2u));
  TrackingVolume* childPtr = child.get();
  childPtr->setNavigationPolicy(
      navPolicyFactory->build(gctx, *childPtr, *logger));
  root->addVolume(std::move(child));

  auto trivialBounds = std::make_shared<const RectangleBounds>(1., 1.);
  auto trivialSurface =
      Surface::makeShared<PlaneSurface>(Transform3::Identity(), trivialBounds);
  trivialSurface->assignGeometryId(
      GeometryIdentifier{}.withVolume(1u).withExtra(11u));
  auto trivialLink =
      std::make_unique<TrivialPortalLink>(trivialSurface, *childPtr);
  root->addPortal(std::make_shared<Portal>(Direction::AlongNormal(),
                                           std::move(trivialLink)));

  auto compositeBounds = std::make_shared<const RectangleBounds>(1., 1.);
  Transform3 compositeTransformA = Transform3::Identity();
  compositeTransformA.pretranslate(Vector3{-1., 0., 0.});
  Transform3 compositeTransformB = Transform3::Identity();
  compositeTransformB.pretranslate(Vector3{1., 0., 0.});
  auto compositeSurfaceA =
      Surface::makeShared<PlaneSurface>(compositeTransformA, compositeBounds);
  auto compositeSurfaceB =
      Surface::makeShared<PlaneSurface>(compositeTransformB, compositeBounds);
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
  auto gridLink = GridPortalLink::make(gridSurface, AxisDirection::AxisX,
                                       Axis{AxisBound, -2., 2., 2});
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
  nlohmann::json encoded = converter.trackingVolumeToJson(gctx, *root);
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

  auto decodedGeometry = converter.fromJson(gctx, encodedFromFile);
  BOOST_REQUIRE(decodedGeometry != nullptr);
  BOOST_REQUIRE(decodedGeometry->highestTrackingVolume() != nullptr);
  BOOST_CHECK_EQUAL(decodedGeometry->highestTrackingVolume()->volumeName(),
                    "root");

  const auto* htvSource = root.get();
  const auto* htvDecoded = decodedGeometry->highestTrackingVolume();
  BOOST_CHECK(*root == *htvDecoded);
  checkHierarchy(gctx, htvSource->volumes(), htvDecoded->volumes());
}

BOOST_AUTO_TEST_CASE(TrackingGeometryJsonConverterNavigation) {
  GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();
  MagneticFieldContext mctx;

  auto field = std::make_shared<ConstantBField>(Vector3{0, 0, 0_T});
  EigenStepper<> stepper{field};

  BoundTrackParameters start = BoundTrackParameters::createCurvilinear(
      Vector4::Zero(), Vector3::UnitX(), 1. / 1_GeV, std::nullopt,
      ParticleHypothesis::pion());

  // Action list and abort list
  using EndOfWorld = EndOfWorldReached;
  using ReferenceActorList = ActorList<StateCollector<>, EndOfWorld>;
  using PropagatorOptions =
      typename ConstantFieldPropagator::template Options<ReferenceActorList>;

  // Options definition
  PropagatorOptions options(gctx, mctx);

  // Build geometry
  CylindricalTrackingGeometry cylindricalGeometryBuilder(gctx, true);
  auto sourceGeometry = cylindricalGeometryBuilder();

  TrackingGeometryJsonConverter converter;
  nlohmann::json encoded = converter.toJson(gctx, *sourceGeometry);
  // TemporaryDirectory tmpDir{};
  auto jsonPath = "tracking_geometry_roundtrip.json";
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

  auto decodedGeometry = converter.fromJson(gctx, encodedFromFile);

  const auto* htvSource = sourceGeometry->highestTrackingVolume();
  const auto* htvDecoded = decodedGeometry->highestTrackingVolume();
  BOOST_CHECK(*htvSource == *htvDecoded);

  checkHierarchy(gctx, htvSource->volumes(), htvDecoded->volumes());

  // -----------------------------------------------------
  // Compare propagation
  Navigator::Config sourceNavCfg;
  sourceNavCfg.trackingGeometry =
      std::make_shared<TrackingGeometry>(*sourceGeometry);
  sourceNavCfg.resolveSensitive = true;
  sourceNavCfg.resolveMaterial = true;
  sourceNavCfg.resolvePassive = false;
  Navigator sourceNavigator{
      sourceNavCfg,
      Acts::getDefaultLogger("SourceNavigator", Acts::Logging::INFO)};
  ConstantFieldPropagator sourcePropagator(stepper, sourceNavigator);

  auto sourceRes = sourcePropagator.propagate(start, options).value();

  Navigator::Config decodedNavCfg;
  decodedNavCfg.trackingGeometry =
      std::make_shared<TrackingGeometry>(*decodedGeometry);
  decodedNavCfg.resolveSensitive = true;
  decodedNavCfg.resolveMaterial = true;
  decodedNavCfg.resolvePassive = false;
  Navigator decodedNavigator{
      sourceNavCfg,
      Acts::getDefaultLogger("DecodedNavigator", Acts::Logging::INFO)};
  ConstantFieldPropagator decodedPropagator(stepper, decodedNavigator);

  auto decodedRes = decodedPropagator.propagate(start, options).value();

  auto& sourceStates = sourceRes.template get<StateCollector<>::result_type>();
  auto& decodedStates =
      decodedRes.template get<StateCollector<>::result_type>();
  for (std::size_t i = 0; i < sourceStates.navigation.size(); i++) {
    auto& navigationA = sourceStates.navigation.at(i);
    auto& navigationB = decodedStates.navigation.at(i);
    if (!navigationA.navigationBreak) {
      BOOST_CHECK(*navigationA.currentVolume == *navigationB.currentVolume);
      BOOST_CHECK(*navigationA.startVolume == *navigationB.startVolume);
      BOOST_CHECK(*navigationA.currentSurface == *navigationB.currentSurface);
      BOOST_CHECK(*navigationA.startSurface == *navigationB.startSurface);
    } else {
      BOOST_CHECK_EQUAL(navigationA.currentVolume, nullptr);
      BOOST_CHECK_EQUAL(navigationB.currentVolume, nullptr);
    }

    auto& steppingA = sourceStates.stepping.at(i);
    auto& steppingB = decodedStates.stepping.at(i);
    BOOST_CHECK(steppingA.pars == steppingB.pars);
    BOOST_CHECK(steppingA.nSteps == steppingB.nSteps);
    BOOST_CHECK(steppingA.pathAccumulated == steppingB.pathAccumulated);
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(MultiLayerNavigationPolicySuite)

namespace {

auto makeMultiLayerVolume() {
  auto tVolume = std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(20., 20., 5.), "CuboidVolume");

  auto boundsRect = std::make_shared<Acts::RectangleBounds>(2., 2.);
  for (int ix = -1; ix <= 1; ++ix) {
    for (int iy = -1; iy <= 1; ++iy) {
      Transform3 trf = Transform3::Identity();
      trf.translation() = Acts::Vector3(4. * ix, 4. * iy, 0.);
      auto surface =
          Acts::Surface::makeShared<Acts::PlaneSurface>(trf, boundsRect);
      surface->assignIsSensitive(true);
      tVolume->addSurface(surface);
    }
  }
  return tVolume;
}

auto makeMultiLayerPolicy(const GeometryContext& gctx,
                          const TrackingVolume& volume, const Logger& logger) {
  Acts::Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisX(-10., 10.,
                                                                   5);
  Acts::Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisY(-10., 10.,
                                                                   5);
  Acts::Grid gridXY(Acts::Type<std::vector<std::size_t>>, std::move(axisX),
                    std::move(axisY));

  Experimental::MultiLayerNavigationPolicy::IndexedUpdatorType indexedGrid(
      std::move(gridXY), {AxisX, AxisY});

  Experimental::MultiLayerNavigationPolicy::Config cfg;
  cfg.binExpansion = {1u, 1u};
  return Experimental::MultiLayerNavigationPolicy(gctx, volume, logger, cfg,
                                                  std::move(indexedGrid));
}

}  // namespace

BOOST_AUTO_TEST_CASE(MultiLayerNavigationPolicyToJson) {
  auto tContext = GeometryContext::dangerouslyDefaultConstruct();
  auto tLogger =
      Acts::getDefaultLogger("MultiLayerNavigationJsonTest", Logging::INFO);

  auto tVolume = makeMultiLayerVolume();
  auto policy = makeMultiLayerPolicy(tContext, *tVolume, *tLogger);

  TrackingGeometryJsonConverter conv;
  nlohmann::json j = conv.navigationPolicyToJson(policy);

  BOOST_CHECK(j.contains("kind"));
  BOOST_CHECK(j.contains("axes"));
  BOOST_CHECK(j.contains("casts"));
  BOOST_CHECK(j.contains("binExpansion"));
  BOOST_CHECK_EQUAL(j["kind"], "MultiLayerNavigation");
  BOOST_CHECK_EQUAL(j["axes"].size(), 2u);
  BOOST_CHECK_EQUAL(j["casts"].size(), 2u);
  BOOST_CHECK(!j.contains("surfaces"));
}

BOOST_AUTO_TEST_CASE(MultiLayerNavigationPolicyRoundTrip) {
  auto tContext = GeometryContext::dangerouslyDefaultConstruct();
  auto tLogger =
      Acts::getDefaultLogger("MultiLayerNavigationJsonTest", Logging::INFO);

  auto tVolume = makeMultiLayerVolume();
  auto policy = makeMultiLayerPolicy(tContext, *tVolume, *tLogger);

  TrackingGeometryJsonConverter conv;
  nlohmann::json j = conv.navigationPolicyToJson(policy);

  auto policyPtr =
      conv.navigationPolicyFromJson(tContext, j, *tVolume, *tLogger);
  auto& restored =
      dynamic_cast<Experimental::MultiLayerNavigationPolicy&>(*policyPtr);

  const auto& origGrid = policy.indexedGrid().grid;
  const auto& restGrid = restored.indexedGrid().grid;
  BOOST_CHECK_EQUAL(origGrid.axes()[0]->getNBins(),
                    restGrid.axes()[0]->getNBins());
  BOOST_CHECK_EQUAL(origGrid.axes()[1]->getNBins(),
                    restGrid.axes()[1]->getNBins());
  BOOST_CHECK_CLOSE(origGrid.axes()[0]->getMin(), restGrid.axes()[0]->getMin(),
                    1e-6);
  BOOST_CHECK_CLOSE(origGrid.axes()[0]->getMax(), restGrid.axes()[0]->getMax(),
                    1e-6);

  BOOST_CHECK_EQUAL(policy.indexedGrid().casts[0],
                    restored.indexedGrid().casts[0]);
  BOOST_CHECK_EQUAL(policy.indexedGrid().casts[1],
                    restored.indexedGrid().casts[1]);

  BOOST_CHECK(policy.config().binExpansion == restored.config().binExpansion);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
