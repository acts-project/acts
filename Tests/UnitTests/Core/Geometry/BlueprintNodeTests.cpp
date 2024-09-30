// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/RootBlueprintNode.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

#include <random>

using namespace Acts::UnitLiterals;

namespace Acts::Test {

auto logger = Acts::getDefaultLogger("UnitTests", Acts::Logging::DEBUG);

GeometryContext gctx;

BOOST_AUTO_TEST_SUITE(Geometry);

BOOST_AUTO_TEST_SUITE(BlueprintNodeTest);

BOOST_AUTO_TEST_CASE(StaticBlueprintNodeConstruction) {
  RootBlueprintNode::Config cfg;
  cfg.envelope[BinningValue::binZ] = {20_mm, 2_mm};
  RootBlueprintNode root{cfg};

  ActsScalar hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);
  auto cyl = std::make_unique<TrackingVolume>(Transform3::Identity(), cylBounds,
                                              "root");

  const auto* cylPtr = cyl.get();

  auto node = std::make_unique<StaticBlueprintNode>(std::move(cyl));

  BOOST_CHECK_EQUAL(node->name(), "root");

  BOOST_CHECK_EQUAL(&node->build({}), cylPtr);

  ObjVisualization3D vis;
  // Add some children
  ActsScalar z0 = -200_mm;
  for (std::size_t i = 0; i < 10; i++) {
    auto childCyl = std::make_unique<TrackingVolume>(
        Transform3::Identity() *
            Translation3{Vector3{0, 0, z0 + i * 2 * hlZ * 1.2}},
        cylBounds, "child" + std::to_string(i));

    GeometryView3D::drawVolume(vis, *childCyl, gctx);

    node->addChild(std::make_unique<StaticBlueprintNode>(std::move(childCyl)));
  }

  BOOST_CHECK_EQUAL(node->children().size(), 10);

  root.addChild(std::move(node));

  std::ofstream ofs{"static.obj"};
  vis.write(ofs);

  auto tGeometry = root.construct({}, gctx, *logger);

  BOOST_REQUIRE(tGeometry);

  BOOST_CHECK_EQUAL(tGeometry->highestTrackingVolume()->volumes().size(), 1);
  std::size_t nVolumes = 0;
  tGeometry->visitVolumes(
      [&](const TrackingVolume* /*volume*/) { nVolumes++; });

  BOOST_CHECK_EQUAL(nVolumes, 12);
}

BOOST_AUTO_TEST_CASE(CylinderContainerNode) {
  ActsScalar hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);

  RootBlueprintNode::Config cfg;
  cfg.envelope[BinningValue::binZ] = {20_mm, 20_mm};
  cfg.envelope[BinningValue::binR] = {0_mm, 20_mm};
  auto root = std::make_unique<RootBlueprintNode>(cfg);

  auto& cyl = root->addCylinderContainer("Container", BinningValue::binZ);

  ActsScalar z0 = -200_mm;
  for (std::size_t i = 0; i < 10; i++) {
    auto childCyl = std::make_unique<TrackingVolume>(
        Transform3::Identity() *
            Translation3{Vector3{0, 0, z0 + i * 2 * hlZ * 1.2}},
        cylBounds, "child" + std::to_string(i));
    cyl.addStaticVolume(std::move(childCyl));
  }

  root->construct({}, gctx, *logger);
}

void pseudoNavigation(const TrackingGeometry& trackingGeometry,
                      Vector3 position, const Vector3& direction,
                      std::ostream& csv, std::size_t run,
                      std::size_t substepsPerCm = 1) {
  std::mt19937 rng{static_cast<unsigned int>(run)};
  std::uniform_real_distribution<> dist{0.01, 0.99};

  const auto* volume = trackingGeometry.lowestTrackingVolume(gctx, position);
  BOOST_REQUIRE_NE(volume, nullptr);
  std::cout << volume->volumeName() << std::endl;

  NavigationStream main;
  const TrackingVolume* currentVolume = volume;

  csv << run << "," << position[0] << "," << position[1] << "," << position[2];
  csv << "," << volume->geometryId().volume();
  csv << "," << volume->geometryId().boundary();
  csv << std::endl;

  std::cout << "start pseudo navigation" << std::endl;

  for (std::size_t i = 0; i < 100; i++) {
    main = NavigationStream{};

    currentVolume->updateNavigationState(
        {.main = main, .position = position, .direction = direction});

    std::cout << main.candidates().size() << " candidates" << std::endl;

    for (const auto& candidate : main.candidates()) {
      std::cout << " -> " << candidate.surface().geometryId() << std::endl;
      std::cout << "    " << candidate.surface().toStream(gctx) << std::endl;
    }

    std::cout << "initializing candidates" << std::endl;
    main.initialize(gctx, {position, direction}, BoundaryTolerance::None());

    std::cout << main.candidates().size() << " candidates remaining"
              << std::endl;

    for (const auto& candidate : main.candidates()) {
      std::cout << " -> " << candidate.surface().geometryId() << std::endl;
      std::cout << "    " << candidate.surface().toStream(gctx) << std::endl;
    }

    if (main.currentCandidate().surface().isOnSurface(gctx, position,
                                                      direction)) {
      std::cout << "Already on portal at initialization, skipping candidate"
                << std::endl;

      auto id = main.currentCandidate().surface().geometryId();
      csv << run << "," << position[0] << "," << position[1] << ","
          << position[2];
      csv << "," << id.volume();
      csv << "," << id.boundary();
      csv << std::endl;
      if (!main.switchToNextCandidate()) {
        std::cout << "candidates exhausted unexpectedly" << std::endl;
        break;
      }
    }

    const auto& candidate = main.currentCandidate();
    std::cout << candidate.portal << std::endl;
    std::cout << candidate.intersection.position().transpose() << std::endl;

    BOOST_REQUIRE_NE(candidate.portal, nullptr);

    std::cout << "on portal: " << candidate.portal->surface().geometryId()
              << std::endl;

    std::cout << "moving to position: " << position.transpose()
              << " (r=" << VectorHelpers::perp(position) << ")\n";
    Vector3 delta = candidate.intersection.position() - position;

    std::size_t substeps =
        std::max(1l, std::lround(delta.norm() / 10_cm * substepsPerCm));

    for (std::size_t j = 0; j < substeps; j++) {
      // position += delta / (substeps + 1);
      Vector3 subpos = position + dist(rng) * delta;
      csv << run << "," << subpos[0] << "," << subpos[1] << "," << subpos[2];
      csv << "," << currentVolume->geometryId().volume();
      csv << "," << currentVolume->geometryId().boundary();
      csv << std::endl;
    }

    position = candidate.intersection.position();

    std::cout << "                 -> " << position.transpose()
              << " (r=" << VectorHelpers::perp(position) << ")" << std::endl;

    currentVolume =
        candidate.portal->resolveVolume(gctx, position, direction).value();

    if (currentVolume == nullptr) {
      std::cout << "switched to nullptr" << std::endl;
      break;
    }

    std::cout << "switched to " << currentVolume->volumeName() << std::endl;

    std::cout << "-----" << std::endl;
  }
}

void portalSamples(const TrackingGeometry& trackingGeometry, Vector3 position,
                   const Vector3& direction, std::ostream& csv,
                   std::size_t run) {
  std::set<const Surface*> visitedSurfaces;

  trackingGeometry.visitVolumes([&](const TrackingVolume* volume) {
    for (const auto& portal : volume->portals()) {
      if (visitedSurfaces.contains(&portal.surface())) {
        continue;
      }
      visitedSurfaces.insert(&portal.surface());

      auto multiIntersection = portal.surface().intersect(
          gctx, position, direction, BoundaryTolerance::None());

      for (const auto& intersection : multiIntersection.split()) {
        if (intersection.isValid()) {
          Vector3 newPosition = intersection.position();
          csv << run << "," << position[0] << "," << position[1] << ","
              << position[2];
          csv << "," << volume->geometryId().volume();
          csv << "," << volume->geometryId().boundary();
          csv << std::endl;
          csv << run << "," << newPosition[0] << "," << newPosition[1] << ","
              << newPosition[2];
          csv << "," << portal.surface().geometryId().volume();
          csv << "," << portal.surface().geometryId().boundary();
          csv << std::endl;
          position = newPosition;
        }
      }
    }
  });
}

BOOST_AUTO_TEST_CASE(NodeApiTestContainers) {
  // Transform3 base{AngleAxis3{30_degree, Vector3{1, 0, 0}}};
  Transform3 base{Transform3::Identity()};

  RootBlueprintNode::Config cfg;
  cfg.envelope[BinningValue::binZ] = {20_mm, 20_mm};
  cfg.envelope[BinningValue::binR] = {0_mm, 20_mm};
  auto root = std::make_unique<RootBlueprintNode>(cfg);

  root->addMaterial([&](auto& mat) {
    mat.addCylinderContainer("Detector", BinningValue::binR, [&](auto& det) {
      det.addCylinderContainer("Pixel", BinningValue::binZ, [&](auto& cyl) {
        cyl.setAttachmentStrategy(CylinderVolumeStack::AttachmentStrategy::Gap)
            .setResizeStrategy(CylinderVolumeStack::ResizeStrategy::Gap);

        cyl.addCylinderContainer(
            "PixelNegativeEndcap", BinningValue::binZ, [&](auto& ec) {
              ec.setAttachmentStrategy(
                  CylinderVolumeStack::AttachmentStrategy::Gap);

              ec.addStaticVolume(
                  base * Translation3{Vector3{0, 0, -600_mm}},
                  std::make_shared<CylinderVolumeBounds>(200_mm, 450_mm, 20_mm),
                  "PixelNeg1");

              ec.addStaticVolume(
                  base * Translation3{Vector3{0, 0, -400_mm}},
                  std::make_shared<CylinderVolumeBounds>(200_mm, 800_mm, 20_mm),
                  "PixelNeg2");
            });

        cyl.addCylinderContainer(
            "PixelBarrel", BinningValue::binR, [&](auto& brl) {
              brl.setAttachmentStrategy(
                     CylinderVolumeStack::AttachmentStrategy::Gap)
                  .setResizeStrategy(
                      CylinderVolumeStack::ResizeStrategy::Expand);

              brl.addStaticVolume(
                  base,
                  std::make_shared<CylinderVolumeBounds>(23_mm, 48_mm, 200_mm),
                  "PixelLayer0");

              brl.addStaticVolume(
                  base,
                  std::make_shared<CylinderVolumeBounds>(87_mm, 103_mm, 250_mm),
                  "PixelLayer1");

              brl.addStaticVolume(base,
                                  std::make_shared<CylinderVolumeBounds>(
                                      150_mm, 180_mm, 310_mm),
                                  "PixelLayer2");

              brl.addStaticVolume(base,
                                  std::make_shared<CylinderVolumeBounds>(
                                      250_mm, 400_mm, 310_mm),
                                  "PixelLayer3");
            });

        auto& ec =
            cyl.addCylinderContainer("PixelPosWrapper", BinningValue::binR);
        ec.setResizeStrategy(CylinderVolumeStack::ResizeStrategy::Gap);
        ec.addStaticVolume(std::make_unique<TrackingVolume>(
            base * Translation3{Vector3{0, 0, 600_mm}},
            std::make_shared<CylinderVolumeBounds>(150_mm, 390_mm, 200_mm),
            "PixelPositiveEndcap"));
      });

      det.addStaticVolume(
          base, std::make_shared<CylinderVolumeBounds>(0_mm, 23_mm, 1000_mm),
          "BeamPipe");
    });
  });

  std::ofstream dot{"api_test_container.dot"};
  root->graphViz(dot);

  auto trackingGeometry = root->construct({}, gctx, *logger);

  trackingGeometry->visitVolumes([&](const TrackingVolume* volume) {
    std::cout << volume->volumeName() << std::endl;
    std::cout << " -> id: " << volume->geometryId() << std::endl;
    std::cout << " -> " << volume->portals().size() << " portals" << std::endl;
  });

  ObjVisualization3D vis;

  trackingGeometry->visualize(vis, gctx, {}, {});

  vis.write("api_test_container.obj");

  Vector3 position = Vector3::Zero();
  std::ofstream csv{"api_test_container.csv"};
  csv << "x,y,z,volume,boundary" << std::endl;

  std::mt19937 rnd{42};

  std::uniform_real_distribution<> dist{-1, 1};

  double etaWidth = 5.0;
  double thetaMin = 2 * std::atan(std::exp(-etaWidth));
  double thetaMax = 2 * std::atan(std::exp(etaWidth));
  std::uniform_real_distribution<> thetaDist{thetaMin, thetaMax};

  using namespace Acts::UnitLiterals;

  for (std::size_t i = 0; i < 5000; i++) {
    // double eta = 1.0;
    // double theta = 2 * std::atan(std::exp(-eta));

    double theta = thetaDist(rnd);
    double phi = 2 * M_PI * dist(rnd);

    Vector3 direction;
    direction[0] = std::sin(theta) * std::cos(phi);
    direction[1] = std::sin(theta) * std::sin(phi);
    direction[2] = std::cos(theta);

    std::cout << "start navigation " << i << std::endl;
    std::cout << "dir: " << direction.transpose() << std::endl;
    std::cout << direction.norm() << std::endl;

    pseudoNavigation(*trackingGeometry, position, direction, csv, i, 2);

    // portalSamples(*trackingGeometry, position, direction, csv, i);

    // break;
  }
}

BOOST_AUTO_TEST_CASE(NodeApiTestConfined) {
  // Transform3 base{AngleAxis3{30_degree, Vector3{1, 0, 0}}};
  Transform3 base{Transform3::Identity()};

  RootBlueprintNode::Config cfg;
  cfg.envelope[BinningValue::binZ] = {20_mm, 20_mm};
  cfg.envelope[BinningValue::binR] = {0_mm, 20_mm};
  auto root = std::make_unique<RootBlueprintNode>(cfg);

  root->addCylinderContainer("Detector", BinningValue::binR, [&](auto& det) {
    det.addStaticVolume(
        base, std::make_shared<CylinderVolumeBounds>(50_mm, 400_mm, 1000_mm),
        "PixelWrapper", [&](auto& wrap) {
          ActsScalar rMin = 100_mm;
          ActsScalar rMax = 350_mm;
          ActsScalar hlZ = 100_mm;

          wrap.addStaticVolume(
              base * Translation3{Vector3{0, 0, -600_mm}},
              std::make_shared<CylinderVolumeBounds>(rMin, rMax, hlZ),
              "PixelNeg1");

          wrap.addStaticVolume(
              base * Translation3{Vector3{0, 0, -200_mm}},
              std::make_shared<CylinderVolumeBounds>(rMin, rMax, hlZ),
              "PixelNeg2");

          wrap.addStaticVolume(
              base * Translation3{Vector3{0, 0, 200_mm}},
              std::make_shared<CylinderVolumeBounds>(rMin, rMax, hlZ),
              "PixelPos1");

          wrap.addStaticVolume(
              base * Translation3{Vector3{0, 0, 600_mm}},
              std::make_shared<CylinderVolumeBounds>(rMin, rMax, hlZ),
              "PixelPos2");
        });

    det.addStaticVolume(
        base, std::make_shared<CylinderVolumeBounds>(0_mm, 23_mm, 1000_mm),
        "BeamPipe");
  });

  std::ofstream dot{"api_test_confined.dot"};
  root->graphViz(dot);

  auto trackingGeometry = root->construct({}, gctx, *logger);

  trackingGeometry->visitVolumes([&](const TrackingVolume* volume) {
    std::cout << volume->volumeName() << std::endl;
    std::cout << " -> id: " << volume->geometryId() << std::endl;
    std::cout << " -> " << volume->portals().size() << " portals" << std::endl;
  });

  ObjVisualization3D vis;

  trackingGeometry->visualize(vis, gctx, {}, {});

  vis.write("api_test_confined.obj");

  const auto* wrapper =
      trackingGeometry->findVolume(GeometryIdentifier{}.setVolume(2));
  BOOST_REQUIRE_NE(wrapper, nullptr);

  std::cout << wrapper->volumeName() << std::endl;

  BOOST_CHECK_EQUAL(wrapper->portals().size(), 20);
  BOOST_CHECK_EQUAL(wrapper->volumes().size(), 4);
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();

}  // namespace Acts::Test
