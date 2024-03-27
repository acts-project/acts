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
  auto cyl =
      std::make_unique<TrackingVolume>(Transform3::Identity(), cylBounds);

  const auto* cylPtr = cyl.get();

  StaticBlueprintNode node("root", std::move(cyl));

  BOOST_CHECK_EQUAL(node.name(), "root");
  node.setName("newroot");
  BOOST_CHECK_EQUAL(node.name(), "newroot");

  BOOST_CHECK_EQUAL(&node.build(), cylPtr);

  ObjVisualization3D vis;
  // Add some children
  ActsScalar z0 = -200_mm;
  for (std::size_t i = 0; i < 10; i++) {
    auto childCyl = std::make_unique<TrackingVolume>(
        Transform3::Identity() *
            Translation3{Vector3{0, 0, z0 + i * 2 * hlZ * 1.2}},
        cylBounds);

    GeometryView3D::drawVolume(vis, *childCyl, gctx);

    node.addChild(std::make_unique<StaticBlueprintNode>(
        "child" + std::to_string(i), std::move(childCyl)));
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
      "root", BinningValue::binZ);

  ActsScalar z0 = -200_mm;
  for (std::size_t i = 0; i < 10; i++) {
    auto childCyl = std::make_unique<TrackingVolume>(
        Transform3::Identity() *
            Translation3{Vector3{0, 0, z0 + i * 2 * hlZ * 1.2}},
        cylBounds);
    root->addChild(std::make_unique<StaticBlueprintNode>(
        "child" + std::to_string(i), std::move(childCyl)));
  }

  // auto world = root->build();
  // std::cout << world.volumeBounds() << std::endl;
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();

}  // namespace Acts::Test
