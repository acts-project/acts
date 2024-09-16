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
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

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

  BOOST_CHECK_EQUAL(&node->build(), cylPtr);

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

  auto world = root.construct(gctx, *logger);

  BOOST_REQUIRE(world);

  BOOST_CHECK_EQUAL(world->highestTrackingVolume()->volumes().size(), 10);
}

BOOST_AUTO_TEST_CASE(CylinderContainerNode) {
  ActsScalar hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);

  RootBlueprintNode::Config cfg;
  auto root = std::make_unique<RootBlueprintNode>(cfg);

  ActsScalar z0 = -200_mm;
  for (std::size_t i = 0; i < 10; i++) {
    auto childCyl = std::make_unique<TrackingVolume>(
        Transform3::Identity() *
            Translation3{Vector3{0, 0, z0 + i * 2 * hlZ * 1.2}},
        cylBounds, "child" + std::to_string(i));
    root->addStaticVolume(std::move(childCyl));
  }

  // BOOST_CHECK_THROW(root->connect(gctx), std::runtime_error);

  TrackingVolume dummy{Transform3::Identity(), cylBounds};
  // BOOST_CHECK_THROW(root->finalize(gctx, dummy), std::runtime_error);

  root->construct(gctx, *logger);
}

BOOST_AUTO_TEST_CASE(NodeApiTestContainers) {
  Transform3 base{AngleAxis3{30_degree, Vector3{1, 0, 0}}};

  RootBlueprintNode::Config cfg;
  cfg.envelope[BinningValue::binZ] = {20_mm, 20_mm};
  cfg.envelope[BinningValue::binR] = {0_mm, 20_mm};
  auto root = std::make_unique<RootBlueprintNode>(cfg);

  root->addMaterial([&](auto& mat) {
    mat.addCylinderContainer("Detector", BinningValue::binR, [&](auto& det) {
      det.addCylinderContainer("Pixel", BinningValue::binZ, [&](auto& cyl) {
        cyl.setAttachmentStrategy(CylinderVolumeStack::AttachmentStrategy::Gap)
            .setResizeStrategy(CylinderVolumeStack::ResizeStrategy::Gap);

        // cyl->addStaticVolume(std::make_unique<TrackingVolume>(
        // base * Translation3{Vector3{0, 0, -600_mm}},
        // std::make_shared<CylinderVolumeBounds>(200_mm, 400_mm, 200_mm),
        // "PixelNegativeEndcap"));

        // cyl->addStaticVolume(std::make_unique<TrackingVolume>(
        // base * Translation3{Vector3{0, 0, 600_mm}},
        // std::make_shared<CylinderVolumeBounds>(200_mm, 400_mm, 200_mm),
        // "PixelPositiveEndcap"));

        // cyl->addStaticVolume(std::make_unique<TrackingVolume>(
        // base * Translation3{Vector3{0, 0, 0_mm}},
        // std::make_shared<CylinderVolumeBounds>(100_mm, 600_mm, 200_mm),
        // "PixelBarrel"));

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

  std::ofstream dot{"api_test.dot"};
  root->graphViz(dot);

  auto trackingGeometry = root->construct(gctx, *logger);

  trackingGeometry->visitVolumes([](const TrackingVolume* volume) {
    std::cout << volume->volumeName() << std::endl;
    std::cout << " -> id: " << volume->geometryId() << std::endl;
    std::cout << " -> " << volume->portals().size() << " portals" << std::endl;
  });

  // ObjVisualization3D vis;
  // root->visualize(vis, gctx);
  // vis.write("api_test.obj");
}

BOOST_AUTO_TEST_CASE(NodeApiTestConfined) {
  Transform3 base{AngleAxis3{30_degree, Vector3{1, 0, 0}}};

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

  auto trackingGeometry = root->construct(gctx, *logger);

  ObjVisualization3D vis;

  trackingGeometry->visitVolumes([&](const TrackingVolume* volume) {
    std::cout << volume->volumeName() << std::endl;
    std::cout << " -> id: " << volume->geometryId() << std::endl;
    std::cout << " -> " << volume->portals().size() << " portals" << std::endl;

    GeometryView3D::drawVolume(vis, *volume, gctx, Transform3::Identity(),
                               {.color = {123, 123, 123}});
  });

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
