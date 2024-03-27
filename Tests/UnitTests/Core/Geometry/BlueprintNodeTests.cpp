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
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

using namespace Acts::UnitLiterals;

namespace Acts::Test {

auto logger = Acts::getDefaultLogger("UnitTests", Acts::Logging::VERBOSE);

GeometryContext gctx;

BOOST_AUTO_TEST_SUITE(Geometry);

BOOST_AUTO_TEST_SUITE(BlueprintNodeTest);

BOOST_AUTO_TEST_CASE(StaticBlueprintNodeConstruction) {
  ActsScalar hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);
  auto cyl = std::make_unique<TrackingVolume>(Transform3::Identity(), cylBounds,
                                              "root");

  const auto* cylPtr = cyl.get();

  StaticBlueprintNode node(std::move(cyl));

  BOOST_CHECK_EQUAL(node.name(), "root");

  BOOST_CHECK_EQUAL(&node.build(), cylPtr);

  ObjVisualization3D vis;
  // Add some children
  ActsScalar z0 = -200_mm;
  for (std::size_t i = 0; i < 10; i++) {
    auto childCyl = std::make_unique<TrackingVolume>(
        Transform3::Identity() *
            Translation3{Vector3{0, 0, z0 + i * 2 * hlZ * 1.2}},
        cylBounds, "child" + std::to_string(i));

    GeometryView3D::drawVolume(vis, *childCyl, gctx);

    node.addChild(std::make_unique<StaticBlueprintNode>(std::move(childCyl)));
  }

  BOOST_CHECK_EQUAL(node.children().size(), 10);

  std::ofstream ofs{"static.obj"};
  vis.write(ofs);

  node.build();
  node.connect();

  auto world = node.releaseVolume();

  BOOST_CHECK_EQUAL(world->volumes().size(), 10);
}

BOOST_AUTO_TEST_CASE(CylinderContainerNode) {
  ActsScalar hlZ = 30_mm;
  auto cylBounds = std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, hlZ);

  auto root = std::make_unique<CylinderContainerBlueprintNode>(
      "root", BinningValue::binZ, CylinderVolumeStack::AttachmentStrategy::Gap,
      CylinderVolumeStack::ResizeStrategy::Gap);

  ActsScalar z0 = -200_mm;
  for (std::size_t i = 0; i < 10; i++) {
    auto childCyl = std::make_unique<TrackingVolume>(
        Transform3::Identity() *
            Translation3{Vector3{0, 0, z0 + i * 2 * hlZ * 1.2}},
        cylBounds, "child" + std::to_string(i));
    root->addStaticVolume(std::move(childCyl));
  }

  TrackingVolume dummy{Transform3::Identity(), cylBounds};
  BOOST_CHECK_THROW(root->connect(dummy), std::runtime_error);

  ObjVisualization3D vis;
  // Can't visualize before having called build
  BOOST_CHECK_THROW(root->visualize(vis, gctx), std::runtime_error);

  auto world = root->build(*logger);

  root->visualize(vis, gctx);

  TrackingVolume top{world.transform(), world.volumeBoundsPtr()};

  root->connect(top, *logger);

  vis.write("container.obj");
}

BOOST_AUTO_TEST_CASE(NodeApiTest) {
  auto root = std::make_unique<CylinderContainerBlueprintNode>(
      "root", BinningValue::binZ, CylinderVolumeStack::AttachmentStrategy::Gap,
      CylinderVolumeStack::ResizeStrategy::Gap);

  root->addCylinderContainer("Pixel", binZ, [&](auto cyl) {
    cyl->setAttachmentStrategy(CylinderVolumeStack::AttachmentStrategy::Gap)
        .setResizeStrategy(CylinderVolumeStack::ResizeStrategy::Gap);

    // cyl->addStaticVolume(std::make_unique<TrackingVolume>(
    // Transform3::Identity() * Translation3{Vector3{0, 0, -600_mm}},
    // std::make_shared<CylinderVolumeBounds>(200_mm, 400_mm, 200_mm),
    // "PixelNegativeEndcap"));

    cyl->addCylinderContainer("PixelNegativeEndcap", binZ, [](auto ec) {
      ec->setAttachmentStrategy(CylinderVolumeStack::AttachmentStrategy::Gap);

      ec->addStaticVolume(std::make_unique<TrackingVolume>(
          Transform3::Identity() * Translation3{Vector3{0, 0, -600_mm}},
          std::make_shared<CylinderVolumeBounds>(200_mm, 450_mm, 20_mm),
          "PixelNeg1"));

      ec->addStaticVolume(std::make_unique<TrackingVolume>(
          Transform3::Identity() * Translation3{Vector3{0, 0, -400_mm}},
          std::make_shared<CylinderVolumeBounds>(200_mm, 800_mm, 20_mm),
          "PixelNeg1"));

      return ec;
    });

    // cyl->addStaticVolume(std::make_unique<TrackingVolume>(
    // Transform3::Identity() * Translation3{Vector3{0, 0, 0_mm}},
    // std::make_shared<CylinderVolumeBounds>(100_mm, 600_mm, 200_mm),
    // "PixelBarrel"));

    cyl->addCylinderContainer("PixelBarrel", binR, [](auto brl) {
      brl->setAttachmentStrategy(CylinderVolumeStack::AttachmentStrategy::Gap)
          .setResizeStrategy(CylinderVolumeStack::ResizeStrategy::Expand);

      brl->addStaticVolume(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(23_mm, 48_mm, 200_mm),
          "PixelLayer0"));

      brl->addStaticVolume(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(87_mm, 103_mm, 250_mm),
          "PixelLayer1"));

      brl->addStaticVolume(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(150_mm, 180_mm, 310_mm),
          "PixelLayer2"));

      brl->addStaticVolume(std::make_unique<TrackingVolume>(
          Transform3::Identity(),
          std::make_shared<CylinderVolumeBounds>(250_mm, 400_mm, 310_mm),
          "PixelLayer2"));

      return brl;
    });

    cyl->addCylinderContainer("PixelPoxWrapper", binR, [](auto ec) {
      ec->setResizeStrategy(CylinderVolumeStack::ResizeStrategy::Gap);
      ec->addStaticVolume(std::make_unique<TrackingVolume>(
          Transform3::Identity() * Translation3{Vector3{0, 0, 600_mm}},
          std::make_shared<CylinderVolumeBounds>(150_mm, 390_mm, 200_mm),
          "PixelPositiveEndcap"));
      return ec;
    });

    return cyl;
  });

  auto world = root->build(*logger);
  ObjVisualization3D vis;
  root->visualize(vis, gctx);
  vis.write("api_test.obj");
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();

}  // namespace Acts::Test
