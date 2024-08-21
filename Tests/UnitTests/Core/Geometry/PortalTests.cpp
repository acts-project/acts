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
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

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
  auto vol2 = makeDummyVolume();

  auto cyl1 = Surface::makeShared<CylinderSurface>(
      Transform3{Translation3{Vector3::UnitZ() * -50_mm}}, 50_mm, 100_mm);

  auto cyl2 = Surface::makeShared<CylinderSurface>(
      Transform3{Translation3{Vector3::UnitZ() * 50_mm}}, 50_mm, 100_mm);

  auto portal1 = std::make_shared<Portal>(
      Direction::AlongNormal, std::make_unique<TrivialPortalLink>(cyl1, *vol1));

  BOOST_CHECK_EQUAL(portal1->resolveVolume(gctx, Vector3{50_mm, 0_mm, -50_mm},
                                           Vector3::UnitX()),
                    vol1.get());

  BOOST_CHECK_EQUAL(portal1->resolveVolume(gctx, Vector3{50_mm, 0_mm, -50_mm},
                                           -Vector3::UnitX()),
                    nullptr);

  auto portal2 = std::make_shared<Portal>(Direction::AlongNormal, cyl2, *vol2);
}
BOOST_AUTO_TEST_CASE(Disc) {}

BOOST_AUTO_TEST_SUITE_END()  // Merging

BOOST_AUTO_TEST_SUITE(Fusing)
BOOST_AUTO_TEST_SUITE_END()  // Fusing

BOOST_AUTO_TEST_CASE(ConstructionFromVolume) {
  // - Cylinder
}

BOOST_AUTO_TEST_SUITE_END()  // Portals

BOOST_AUTO_TEST_SUITE_END()  // Geometry

}  // namespace Acts::Test
