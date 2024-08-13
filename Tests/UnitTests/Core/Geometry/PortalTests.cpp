// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/tools/context.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceMergingException.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

using namespace Acts::UnitLiterals;

namespace Acts::Test {

auto logger = Acts::getDefaultLogger("UnitTests", Acts::Logging::VERBOSE);

struct Fixture {
  Logging::Level m_level;
  Fixture() {
    m_level = Acts::Logging::getFailureThreshold();
    Acts::Logging::setFailureThreshold(Acts::Logging::FATAL);
  }

  ~Fixture() { Acts::Logging::setFailureThreshold(m_level); }
};

std::shared_ptr<TrackingVolume> makeDummyVolume() {
  return std::make_shared<TrackingVolume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(30_mm, 40_mm, 100_mm));
}

GeometryContext gctx;

BOOST_FIXTURE_TEST_SUITE(Geometry, Fixture)

BOOST_AUTO_TEST_SUITE(Portals)
BOOST_AUTO_TEST_SUITE(Merging)

BOOST_AUTO_TEST_CASE(Cylinder) {
  auto vol1 = makeDummyVolume();
  vol1->setVolumeName("vol1");
  auto vol2 = makeDummyVolume();
  vol2->setVolumeName("vol2");

  auto cyl1 = Surface::makeShared<CylinderSurface>(
      Transform3{Translation3{Vector3::UnitZ() * -100_mm}}, 50_mm, 100_mm);

  auto cyl2 = Surface::makeShared<CylinderSurface>(
      Transform3{Translation3{Vector3::UnitZ() * 100_mm}}, 50_mm, 100_mm);

  auto portal1 = std::make_shared<Portal>(
      Direction::AlongNormal, std::make_unique<TrivialPortalLink>(cyl1, *vol1));
  BOOST_REQUIRE(portal1);

  BOOST_CHECK_EQUAL(portal1->resolveVolume(gctx, Vector3{50_mm, 0_mm, -50_mm},
                                           Vector3::UnitX()),
                    vol1.get());

  BOOST_CHECK_EQUAL(portal1->resolveVolume(gctx, Vector3{50_mm, 0_mm, -50_mm},
                                           -Vector3::UnitX()),
                    nullptr);

  auto portal2 = std::make_shared<Portal>(Direction::AlongNormal, cyl2, *vol2);
  BOOST_REQUIRE(portal2);

  BOOST_CHECK_EQUAL(portal2->resolveVolume(gctx, Vector3{50_mm, 0_mm, -50_mm},
                                           -Vector3::UnitX()),
                    nullptr);

  BOOST_CHECK_EQUAL(portal2->resolveVolume(gctx, Vector3{50_mm, 0_mm, -50_mm},
                                           Vector3::UnitX()),
                    vol2.get());

  auto portal3 = std::make_shared<Portal>(
      std::make_unique<TrivialPortalLink>(cyl2, *vol2), nullptr);

  BOOST_CHECK_NE(portal3->getLink(Direction::AlongNormal), nullptr);
  BOOST_CHECK_EQUAL(portal3->getLink(Direction::OppositeNormal), nullptr);

  auto portal4 = std::make_shared<Portal>(
      nullptr, std::make_unique<TrivialPortalLink>(cyl2, *vol2));

  BOOST_CHECK_EQUAL(portal4->getLink(Direction::AlongNormal), nullptr);
  BOOST_CHECK_NE(portal4->getLink(Direction::OppositeNormal), nullptr);

  // Not mergeable because 1 has portal along but 4 has portal oppsite
  //         ^
  //         |
  //  portal1|              portal2
  // +-------+-------+  +  +---------------+
  // |               |     |               |
  // +---------------+     +-------+-------+
  //                               |
  //                               |
  //                               v
  BOOST_CHECK_THROW(
      Portal::merge(portal1, portal4, BinningValue::binZ, *logger),
      PortalMergingException);

  BOOST_CHECK_EQUAL(portal2->resolveVolume(gctx, Vector3{50_mm, 0_mm, 50_mm},
                                           Vector3::UnitX()),
                    vol2.get());

  BOOST_CHECK_EQUAL(portal2->resolveVolume(gctx, Vector3{50_mm, 0_mm, 50_mm},
                                           -Vector3::UnitX()),
                    nullptr);

  BOOST_CHECK_THROW(
      Portal::merge(portal1, portal2, BinningValue::binRPhi, *logger),
      SurfaceMergingException);

  //         ^                     ^
  //         |                     |
  //  portal1|              portal2|
  // +-------+-------+  +  +-------+-------+
  // |               |     |               |
  // +---------------+     +---------------+

  auto merged12 = Portal::merge(portal1, portal2, BinningValue::binZ, *logger);
  BOOST_REQUIRE(merged12);
  BOOST_CHECK_NE(merged12->getLink(Direction::AlongNormal), nullptr);
  BOOST_CHECK_EQUAL(merged12->getLink(Direction::OppositeNormal), nullptr);

  auto grid12 = dynamic_cast<const GridPortalLink*>(
      merged12->getLink(Direction::AlongNormal));
  grid12->printContents(std::cout);

  BOOST_CHECK_EQUAL(merged12->resolveVolume(gctx, Vector3{50_mm, 0_mm, -50_mm},
                                            Vector3::UnitX()),
                    vol1.get());

  BOOST_CHECK_EQUAL(merged12->resolveVolume(gctx, Vector3{50_mm, 0_mm, 50_mm},
                                            Vector3::UnitX()),
                    vol2.get());

  BOOST_CHECK_EQUAL(merged12->resolveVolume(gctx, Vector3{50_mm, 0_mm, -50_mm},
                                            -Vector3::UnitX()),
                    nullptr);

  BOOST_CHECK_EQUAL(merged12->resolveVolume(gctx, Vector3{50_mm, 0_mm, 50_mm},
                                            -Vector3::UnitX()),
                    nullptr);

  // Can't merge with self because surfaces are not mergeable
  BOOST_CHECK_THROW(
      Portal::merge(portal1, portal1, BinningValue::binZ, *logger),
      AssertionFailureException);
}

BOOST_AUTO_TEST_CASE(Disc) {
  auto vol1 = makeDummyVolume();
  vol1->setVolumeName("vol1");
  auto vol2 = makeDummyVolume();
  vol2->setVolumeName("vol2");
  auto vol3 = makeDummyVolume();
  vol3->setVolumeName("vol3");
  auto vol4 = makeDummyVolume();
  vol4->setVolumeName("vol4");

  auto disc1 = Surface::makeShared<DiscSurface>(
      Transform3::Identity(), std::make_shared<RadialBounds>(50_mm, 100_mm));

  auto disc2 = Surface::makeShared<DiscSurface>(
      Transform3::Identity(), std::make_shared<RadialBounds>(100_mm, 150_mm));

  auto portal1 = std::make_shared<Portal>(
      std::make_unique<TrivialPortalLink>(disc1, *vol1),
      std::make_unique<TrivialPortalLink>(disc1, *vol2));

  auto portal2 = std::make_shared<Portal>(
      std::make_unique<TrivialPortalLink>(disc2, *vol3),
      std::make_unique<TrivialPortalLink>(disc2, *vol4));

  BOOST_CHECK_EQUAL(portal1->resolveVolume(gctx, Vector3{55_mm, 0_mm, 0_mm},
                                           Vector3::UnitZ()),
                    vol1.get());
  BOOST_CHECK_EQUAL(portal1->resolveVolume(gctx, Vector3{55_mm, 0_mm, 0_mm},
                                           -Vector3::UnitZ()),
                    vol2.get());

  BOOST_CHECK_EQUAL(portal2->resolveVolume(gctx, Vector3{105_mm, 0_mm, 0_mm},
                                           Vector3::UnitZ()),
                    vol3.get());
  BOOST_CHECK_EQUAL(portal2->resolveVolume(gctx, Vector3{105_mm, 0_mm, 0_mm},
                                           -Vector3::UnitZ()),
                    vol4.get());

  BOOST_CHECK_THROW(
      Portal::merge(portal1, portal2, BinningValue::binZ, *logger),
      AssertionFailureException);

  BOOST_CHECK_THROW(
      Portal::merge(portal1, portal2, BinningValue::binPhi, *logger),
      SurfaceMergingException);

  //         ^                     ^
  //         |                     |
  //  portal1|              portal2|
  // +-------+-------+     +-------+-------+
  // |               |  +  |               |
  // +-------+-------+     +-------+-------+
  //         |                     |
  //         |                     |
  //         v                     v
  auto merged12 = Portal::merge(portal1, portal2, BinningValue::binR, *logger);
  BOOST_REQUIRE(merged12);

  BOOST_CHECK_EQUAL(merged12->resolveVolume(gctx, Vector3{55_mm, 0_mm, 0_mm},
                                            Vector3::UnitZ()),
                    vol1.get());
  BOOST_CHECK_EQUAL(merged12->resolveVolume(gctx, Vector3{55_mm, 0_mm, 0_mm},
                                            -Vector3::UnitZ()),
                    vol2.get());

  BOOST_CHECK_EQUAL(merged12->resolveVolume(gctx, Vector3{105_mm, 0_mm, 0_mm},
                                            Vector3::UnitZ()),
                    vol3.get());
  BOOST_CHECK_EQUAL(merged12->resolveVolume(gctx, Vector3{105_mm, 0_mm, 0_mm},
                                            -Vector3::UnitZ()),
                    vol4.get());
}

BOOST_AUTO_TEST_SUITE_END()  // Merging

BOOST_AUTO_TEST_SUITE(Fusing)

// @TODO: Test fusing portals with different surfaces fails
BOOST_AUTO_TEST_CASE(Separated) {
  auto vol1 = makeDummyVolume();
  vol1->setVolumeName("vol1");
  auto vol2 = makeDummyVolume();
  vol2->setVolumeName("vol2");

  auto cyl1 = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                   50_mm, 100_mm);

  auto cyl2 = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                   60_mm, 100_mm);

  auto portal1 = std::make_shared<Portal>(
      Direction::OppositeNormal,
      std::make_unique<TrivialPortalLink>(cyl1, *vol1));
  auto portal2 = std::make_shared<Portal>(
      Direction::AlongNormal, std::make_unique<TrivialPortalLink>(cyl2, *vol2));

  BOOST_CHECK_THROW(Portal::fuse(portal1, portal2, *logger),
                    PortalFusingException);

  auto disc1 = Surface::makeShared<DiscSurface>(
      Transform3::Identity(), std::make_shared<RadialBounds>(50_mm, 100_mm));

  auto disc2 = Surface::makeShared<DiscSurface>(
      Transform3{Translation3{Vector3{0, 0, 5_mm}}},
      std::make_shared<RadialBounds>(50_mm, 100_mm));

  portal1 = std::make_shared<Portal>(
      Direction::OppositeNormal,
      std::make_unique<TrivialPortalLink>(disc1, *vol1));
  portal2 = std::make_shared<Portal>(
      Direction::AlongNormal,
      std::make_unique<TrivialPortalLink>(disc2, *vol2));

  BOOST_CHECK_THROW(Portal::fuse(portal1, portal2, *logger),
                    PortalFusingException);
}

BOOST_AUTO_TEST_SUITE_END()  // Fusing

BOOST_AUTO_TEST_CASE(ConstructionFromVolume) {
  // - Cylinder
}

BOOST_AUTO_TEST_CASE(Construction) {
  // @TODO: Displaced surfaces fail
  // @TODO: Surface unification to same instance
}

BOOST_AUTO_TEST_SUITE_END()  // Portals

BOOST_AUTO_TEST_SUITE_END()  // Geometry

}  // namespace Acts::Test
