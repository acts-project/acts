// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/MultiWireVolumeBuilder.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrapezoidPortalShell.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/TransformHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

#include <cmath>
#include <memory>
#include <numbers>
#include <vector>

using namespace Acts;
using namespace Acts::UnitLiterals;

auto tContext = GeometryContext::dangerouslyDefaultConstruct();

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(NavigationSuite)
ACTS_LOCAL_LOGGER(getDefaultLogger("MultiWireNavigationTests",
                                   Logging::Level::VERBOSE));

// advance pos to the intersection with the current target surface
void step(const GeometryContext& geoCtx, Vector3& pos, const Vector3& dir,
          const NavigationTarget& target) {
  auto isect = target.surface()
                   .intersect(geoCtx, pos, dir, target.boundaryTolerance())
                   .closest();
  pos += isect.pathLength() * dir;
}

// build a grid of staggered straw surfaces for the first test
void generateStrawSurfaces(const TrapezoidVolumeBounds& volBounds,
                           const Transform3& localToGlobal,
                           std::vector<std::shared_ptr<Surface>>& out,
                           ObjVisualization3D& vis) {
  constexpr double strawRadius = 5._cm;
  const double halfX =
      std::max(volBounds.get(TrapezoidVolumeBounds::eHalfLengthXnegY),
               volBounds.get(TrapezoidVolumeBounds::eHalfLengthXposY));
  const double halfY = volBounds.get(TrapezoidVolumeBounds::eHalfLengthY);
  const double halfZ = volBounds.get(TrapezoidVolumeBounds::eHalfLengthZ);

  const auto nLayers = static_cast<std::size_t>(
      std::floor(2 * halfZ / (std::sqrt(3.) * strawRadius)));
  std::size_t nStraws = std::floor(((halfX) / strawRadius));
  std::cout << nLayers << std::endl;
  std::cout << nStraws << std::endl;
  Vector3 ipos = {-halfX + strawRadius, -0., -halfZ + strawRadius};
  Vector3 pos = ipos;
  auto strawBounds = std::make_shared<LineBounds>(strawRadius, halfY - 0.5_mm);

  for (std::size_t i = 0; i < nLayers; i++) {
    for (std::size_t j = 0; j < nStraws; j++) {
      const Vector3 localPos{
          ipos.x() + (2 * j + (i) % 2) * strawRadius,  // x: across the layer
          0,                                           // y: layer-to-layer
          ipos.z() +
              i * std::sqrt(3.) * strawRadius};  // z: along the wire (centered)
      Transform3 tubetrans =
          Transform3(Translation3(localPos)) *
          Transform3(AngleAxis3{90._degree, Vector3::UnitX()});
      Transform3 trans = localToGlobal * tubetrans;
      auto newStraw = out.emplace_back(
          Surface::makeShared<StrawSurface>(trans, strawBounds));
      newStraw->assignIsSensitive(true);
      newStraw->assignGeometryId(
          GeometryIdentifier{}.withLayer(i + 1).withSensitive(j + 1));
      GeometryView3D::drawSurface(vis, *newStraw, tContext,
                                  Transform3::Identity());
    }
  }
}

// Test 1: Candidates initialization test inside a multiwire volume
BOOST_AUTO_TEST_CASE(MultiLayer_NavigationPolicy) {
  std::vector<std::shared_ptr<Surface>> strawSurfaces{};
  ObjVisualization3D visualHelper{};

  auto volBounds =
      std::make_shared<TrapezoidVolumeBounds>(0.925_m, 0.925_m, 1.2_m, 0.14_m);
  const Transform3 volTrans =
      Transform3(Translation3(Vector3(300., -150., 500.))) *  // translation
      Transform3(AngleAxis3(35._degree, Vector3::UnitZ())) *  // rotation
      Transform3(AngleAxis3(20._degree, Vector3::UnitX()));
  std::cout << "here" << std::endl;

  generateStrawSurfaces(*volBounds, volTrans, strawSurfaces, visualHelper);

  MultiWireVolumeBuilder::Config mwCfg;
  mwCfg.name = "MultiWireVolume";
  mwCfg.mlSurfaces = strawSurfaces;
  mwCfg.binning = {{AxisDirection::AxisX, 0u}, {AxisDirection::AxisZ, 0u}};
  mwCfg.shiftDirection = AxisDirection::AxisX;
  mwCfg.bounds = volBounds;
  mwCfg.transform = volTrans;

  std::cout << "surface sgenerated" << std::endl;

  MultiWireVolumeBuilder mwBuilder(mwCfg);
  std::unique_ptr<TrackingVolume> volume = mwBuilder.buildVolume();

  GeometryView3D::drawVolume(visualHelper, *volume, tContext,
                             Transform3::Identity());

  std::cout << "surface sgenerated" << std::endl;

  visualHelper.write("MultiLayerNavigation_test1.obj");

  SingleTrapezoidPortalShell portalShell{tContext, *volume};
  portalShell.applyToVolume();

  BOOST_CHECK(volume->volumes().empty());
  // 18 straws per layer and in total 3 layers
  BOOST_CHECK_EQUAL(volume->surfaces().size(), 54u);
  BOOST_CHECK_EQUAL(volume->portals().size(), 6u);

  NavigationStream main;
  AppendOnlyNavigationStream stream{main};
  Vector3 startPos = {0., 0., -199.};
  Vector3 startDir = {0., 0., 1.};
  NavigationArguments args{startPos, startDir};

  auto navFactory = mwBuilder.createNavigationPolicyFactory(tContext);
  volume->setNavigationPolicy(navFactory->build(tContext, *volume, logger()));

  NavigationPolicyStateManager stateManager;
  volume->navigationPolicy()->createState(tContext, args, stateManager,
                                          logger());
  auto policyState = stateManager.currentState();
  volume->initializeNavigationCandidates(tContext, args, policyState, stream,
                                         logger());
  std::cout << "surface sgenerated" << std::endl;
  BOOST_CHECK_EQUAL(main.candidates().size(), 9u);

  auto it = std::unique(main.candidates().begin(), main.candidates().end(),
                        [](const auto& lhs, const auto& rhs) {
                          return lhs.surface() == rhs.surface();
                        });
  BOOST_CHECK(it == main.candidates().end());

  double angle = std::numbers::pi / 4.;
  startDir = {std::cos(angle), 0., std::sin(angle)};
  args.direction = startDir;

  // clear the candidates and re initialize with new arguments
  main.reset();
  NavigationPolicyStateManager stateManager2;
  volume->navigationPolicy()->createState(tContext, args, stateManager2,
                                          logger());
  auto policyState2 = stateManager2.currentState();
  volume->initializeNavigationCandidates(tContext, args, policyState2, stream,
                                         logger());
  BOOST_CHECK_EQUAL(main.candidates().size(), 9u);
}

// // Test 2: navigate a ring of impact points around each straw's center and
// check if target is reached
BOOST_AUTO_TEST_CASE(MultiLayerNavigation_TargetSurfaces) {
  std::vector<std::shared_ptr<Surface>> straws{};
  ObjVisualization3D vis{};

  auto volBounds =
      std::make_shared<TrapezoidVolumeBounds>(0.925_m, 0.925_m, 1.2_m, 0.14_m);
  const Transform3 volTrans =
      Transform3(Translation3(Vector3(300., -500., 500.))) *  // translation
      Transform3(AngleAxis3(35._degree, Vector3::UnitZ())) *  // rotation
      Transform3(AngleAxis3(20._degree, Vector3::UnitX()));

  generateStrawSurfaces(*volBounds, volTrans, straws, vis);

  Blueprint::Config bpCfg{};
  bpCfg.envelope[AxisDirection::AxisX] = {20_mm, 20_mm};
  bpCfg.envelope[AxisDirection::AxisY] = {20_mm, 20_mm};
  bpCfg.envelope[AxisDirection::AxisZ] = {20_mm, 20_mm};
  Blueprint root{bpCfg};

  auto& container = root.addStaticVolume(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(25._m, 25._m, 25._m), "world");

  // start plane: launch tracks from the origin toward each straw.
  const Transform3 surfaceTrans =
      Transform3(Translation3(Vector3(0., 0., -1000.)));

  auto startSurface = Surface::makeShared<PlaneSurface>(
      surfaceTrans, std::make_shared<RectangleBounds>(10._m, 10._m));

  GeometryView3D::drawSurface(vis, *startSurface, tContext,
                              Transform3::Identity());

  {
    MultiWireVolumeBuilder::Config mwCfg{};
    mwCfg.name = "MultiWireVolume";
    mwCfg.mlSurfaces = straws;
    mwCfg.binning = {{AxisDirection::AxisX, 0u}, {AxisDirection::AxisZ, 0u}};
    mwCfg.shiftDirection = AxisDirection::AxisX;
    mwCfg.bounds = volBounds;
    mwCfg.transform = volTrans;

    MultiWireVolumeBuilder mwBuilder{mwCfg};
    auto volume = mwBuilder.buildVolume();
    container.addStaticVolume(std::move(volume))
        .setNavigationPolicyFactory(
            mwBuilder.createNavigationPolicyFactory(tContext));
  }

  std::shared_ptr<const TrackingGeometry> trkGeo =
      root.construct(BlueprintOptions{}, tContext, logger());
  BOOST_REQUIRE(trkGeo != nullptr);

  const auto* mwVolume = trkGeo->findVolumeByName("MultiWireVolume");
  BOOST_REQUIRE(mwVolume != nullptr);
  BOOST_CHECK(mwVolume->geometryId() != GeometryIdentifier{});

  trkGeo->visitVolumes([&](const TrackingVolume* vol) {
    GeometryView3D::drawVolume(vis, *vol, tContext, Transform3::Identity());
  });

  vis.write("MultiLayerNavTest_trkGeo.obj");

  // navigator configuration
  Navigator::Config navCfg;
  navCfg.trackingGeometry = std::move(trkGeo);
  navCfg.resolveSensitive = true;
  navCfg.resolveMaterial = true;
  navCfg.resolvePassive = false;
  Navigator navigator{navCfg, getDefaultLogger("MWNav", Logging::VERBOSE)};
  Vector3 start = {0., 0., -1000.};
  Vector3 dir = {0., 0., 1.};
  Navigator::Options options{tContext};
  Navigator::State state = navigator.makeState(options);
  NavigationTarget target = navigator.nextTarget(state, start, dir);
  NavigatorInitializeArguments initArgs;
  initArgs.position = start;
  initArgs.direction = dir;
  BOOST_CHECK(navigator.initialize(state, initArgs).ok());
  BOOST_CHECK(!target.isNone());
  // expect to be a boundary
  BOOST_CHECK(target.isPortalTarget());

  // ring of impact points around each straw center
  constexpr std::size_t nSamples = 5;
  constexpr double rFrac = 0.5;
  constexpr double tubeR = 5._cm;
  for (std::size_t s = 0; s < straws.size(); ++s) {
    const auto& straw = straws.at(s);
    const Transform3 l2g = straw->localToGlobalTransform(tContext);

    for (std::size_t k = 0; k < nSamples; ++k) {
      const double phi = 2. * std::numbers::pi * k / nSamples;
      const Vector3 locPos =
          rFrac * tubeR * makeDirectionFromPhiTheta(phi, 90._degree);
      const Vector3 impact = l2g * locPos;
      ACTS_VERBOSE(__LINE__
                   << " - Start to target surface: " << straw->geometryId()
                   << ", impact position " << toString(impact) << ", local: "
                   << toString(locPos) << " phi: " << (phi / 1._degree));
      Vector3 pos = start;
      const Vector3 direction = (impact - pos).normalized();

      Navigator::State st = navigator.makeState(options);
      NavigatorInitializeArguments ia;
      ia.position = pos;
      ia.direction = direction;
      ia.startSurface =
          startSurface.get();  // no start surface; enter from outside
      if (!navigator.initialize(st, ia).ok()) {
        BOOST_CHECK(false);
        continue;
      }

      NavigationTarget tgt = navigator.nextTarget(st, pos, direction);
      ACTS_VERBOSE(__LINE__ << " - The start target is " << tgt);
      while (!tgt.isNone() &&
             tgt.surface().geometryId() != straw->geometryId()) {
        step(tContext, pos, direction, tgt);
        navigator.handleSurfaceReached(st, pos, direction, tgt.surface());
        tgt = navigator.nextTarget(st, pos, direction);
        ACTS_VERBOSE(__LINE__ << " - Continue to next target " << tgt);
      }
      ACTS_VERBOSE(__LINE__ << " - Final target is " << tgt);
      BOOST_CHECK(tgt.isSurfaceTarget());

      BOOST_CHECK_EQUAL(&tgt.surface(), straw.get());
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
