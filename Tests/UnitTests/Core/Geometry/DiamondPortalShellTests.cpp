// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/DiamondPortalShell.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

GeometryContext gctx = GeometryContext::dangerouslyDefaultConstruct();

BOOST_AUTO_TEST_SUITE(GeometrySuite)

BOOST_AUTO_TEST_CASE(ConstructionFromVolume) {
  // check exception throw if not the correct shape

  TrackingVolume fVol(
      Transform3::Identity(),
      std::make_shared<TrapezoidVolumeBounds>(20._cm, 10._cm, 10._cm, 5._cm));

  BOOST_CHECK_THROW(SingleDiamondPortalShell(gctx, fVol),
                    std::invalid_argument);

  // conastruct a convex polygon tracking volume for which we are gonna build
  // the portal shell

  TrackingVolume testVol(Transform3::Identity(),
                         std::make_shared<DiamondVolumeBounds>(
                             20._cm, 25._cm, 15._cm, 15._cm, 20._cm, 12._cm));

  SingleDiamondPortalShell polygShell{gctx, testVol};

  BOOST_CHECK(polygShell.isValid());
  BOOST_CHECK_EQUAL(polygShell.size(), 8);

  // check if the portals are correctly built from the faces
  using enum DiamondVolumeBounds::Face;

  const auto nXY = polygShell.portalPtr(NegativeZFaceXY);
  const auto pXY = polygShell.portalPtr(PositiveZFaceXY);
  const auto nXZ = polygShell.portalPtr(NegativeYFaceZX);
  const auto pXZ = polygShell.portalPtr(PositiveYFaceZX);
  const auto nYZ12 = polygShell.portalPtr(NegativeXFaceYZ12);
  const auto pYZ12 = polygShell.portalPtr(PositiveXFaceYZ12);
  const auto nYZ23 = polygShell.portalPtr(NegativeXFaceYZ23);
  const auto pYZ23 = polygShell.portalPtr(PositiveXFaceYZ23);

  double alphaAngle = std::numbers::pi - std::atan2(15._cm, 5._cm);
  double betaAngle = std::atan2(20._cm, 10._cm);

  BOOST_REQUIRE_NE(nXY, nullptr);
  BOOST_REQUIRE_NE(pXY, nullptr);
  BOOST_REQUIRE_NE(nXZ, nullptr);
  BOOST_REQUIRE_NE(pXZ, nullptr);
  BOOST_REQUIRE_NE(nYZ12, nullptr);
  BOOST_REQUIRE_NE(pYZ12, nullptr);
  BOOST_REQUIRE_NE(nYZ23, nullptr);
  BOOST_REQUIRE_NE(pYZ23, nullptr);

  BOOST_CHECK_EQUAL(
      nXY->resolveVolume(gctx, Vector3{0., 0., -12._cm}, Vector3::UnitZ())
          .value(),
      &testVol);
  BOOST_CHECK_EQUAL(
      pXY->resolveVolume(gctx, Vector3{0., 0., 12._cm}, -Vector3::UnitZ())
          .value(),
      &testVol);
  BOOST_CHECK_EQUAL(
      nXZ->resolveVolume(gctx, Vector3{0., -15._cm, 0.}, Vector3::UnitY())
          .value(),
      &testVol);

  BOOST_CHECK_EQUAL(
      pXZ->resolveVolume(gctx, Vector3{0., 20._cm, 0.}, -Vector3::UnitY())
          .value(),
      &testVol);

  Vector3 normalVec =
      Vector3{std::cos(-std::numbers::pi / 2. + alphaAngle),
              std::sin(-std::numbers::pi / 2. + alphaAngle), 0.};

  BOOST_CHECK_EQUAL(
      nYZ12->resolveVolume(gctx, Vector3{-22.5_cm, -7.5_cm, 0.}, normalVec)
          .value(),
      &testVol);

  BOOST_CHECK_EQUAL(
      pYZ12->resolveVolume(gctx, Vector3{22.5_cm, -7.5_cm, 0.}, -normalVec)
          .value(),
      &testVol);

  normalVec = {std::cos(std::numbers::pi / 2. - betaAngle),
               std::sin(std::numbers::pi / 2. - betaAngle), 0.};
  BOOST_CHECK_EQUAL(
      pYZ23->resolveVolume(gctx, Vector3{20._cm, 10._cm, 0}, -normalVec)
          .value(),
      &testVol);

  BOOST_CHECK_EQUAL(
      nYZ23->resolveVolume(gctx, Vector3{-20._cm, 10._cm, 0.}, normalVec)
          .value(),
      &testVol);
}

BOOST_AUTO_TEST_CASE(PortalAssignment) {
  using enum DiamondVolumeBounds::Face;

  TrackingVolume polygVol(Transform3::Identity(),
                          std::make_shared<DiamondVolumeBounds>(
                              20._cm, 25._cm, 15._cm, 15._cm, 20._cm, 12._cm));

  SingleDiamondPortalShell polygShell{gctx, polygVol};

  // get the portal faces
  const auto nXY = polygShell.portalPtr(NegativeZFaceXY);
  const auto pXY = polygShell.portalPtr(PositiveZFaceXY);

  BOOST_REQUIRE_NE(nXY, nullptr);
  BOOST_REQUIRE_NE(pXY, nullptr);

  const auto nYZ12 = polygShell.portalPtr(NegativeXFaceYZ12);
  const auto pYZ12 = polygShell.portalPtr(PositiveXFaceYZ12);

  const auto nYZ23 = polygShell.portalPtr(NegativeXFaceYZ23);
  const auto pYZ23 = polygShell.portalPtr(PositiveXFaceYZ23);

  const auto nZX = polygShell.portalPtr(NegativeYFaceZX);
  const auto pZX = polygShell.portalPtr(PositiveYFaceZX);

  BOOST_REQUIRE_NE(nYZ23, nullptr);
  BOOST_REQUIRE_NE(pYZ23, nullptr);

  BOOST_REQUIRE_NE(nZX, nullptr);
  BOOST_REQUIRE_NE(pZX, nullptr);

  BOOST_REQUIRE_NE(nYZ12, nullptr);
  BOOST_REQUIRE_NE(pYZ12, nullptr);

  auto* pZXLink = dynamic_cast<const TrivialPortalLink*>(
      pZX->getLink(Direction::OppositeNormal()));
  auto grid = pZXLink->makeGrid(AxisDirection::AxisY);

  BOOST_REQUIRE_NE(pZXLink, nullptr);

  auto newPortalpZX =
      std::make_shared<Portal>(Direction::OppositeNormal(), std::move(grid));
  polygShell.setPortal(newPortalpZX, PositiveYFaceZX);

  BOOST_CHECK_EQUAL(polygShell.portalPtr(PositiveYFaceZX), newPortalpZX);

  // //check also the other portals
  BOOST_CHECK_EQUAL(polygShell.portalPtr(NegativeZFaceXY), nXY);
  BOOST_CHECK_EQUAL(polygShell.portalPtr(PositiveZFaceXY), pXY);
  BOOST_CHECK_EQUAL(polygShell.portalPtr(NegativeXFaceYZ12), nYZ12);
  BOOST_CHECK_EQUAL(polygShell.portalPtr(PositiveXFaceYZ12), pYZ12);
  BOOST_CHECK_EQUAL(polygShell.portalPtr(NegativeXFaceYZ23), nYZ23);
  BOOST_CHECK_EQUAL(polygShell.portalPtr(PositiveXFaceYZ23), pYZ23);
  BOOST_CHECK_EQUAL(polygShell.portalPtr(NegativeYFaceZX), nZX);
}

BOOST_AUTO_TEST_CASE(Fill) {
  using enum DiamondVolumeBounds::Face;

  TrackingVolume testVol(Transform3::Identity(),
                         std::make_shared<DiamondVolumeBounds>(
                             20._cm, 25._cm, 15._cm, 15._cm, 20._cm, 12._cm));
  SingleDiamondPortalShell polygShell{gctx, testVol};

  // without filling the protal shell from a volume the portal link to this
  // direction shouldn't exist - but only the other direction
  BOOST_CHECK_EQUAL(polygShell.portalPtr(NegativeZFaceXY)
                        ->getLink(Direction::OppositeNormal()),
                    nullptr);
  BOOST_CHECK(polygShell.portalPtr(NegativeZFaceXY)
                  ->getLink(Direction::AlongNormal()) != nullptr);

  // create a volume to link to the portal to the other side
  TrackingVolume testVol2(
      Transform3::Identity() * Translation3(Vector3::UnitZ() * 24_cm),
      std::make_shared<DiamondVolumeBounds>(20._cm, 25._cm, 15._cm, 15._cm,
                                            20._cm, 12._cm));
  polygShell.fill(testVol2);

  BOOST_CHECK_NE(polygShell.portalPtr(NegativeZFaceXY)
                     ->getLink(Direction::OppositeNormal()),
                 nullptr);
}

BOOST_AUTO_TEST_CASE(ApplyToVolume) {
  // test the volume's portals
  TrackingVolume testVol(Transform3::Identity(),
                         std::make_shared<DiamondVolumeBounds>(
                             20._cm, 25._cm, 15._cm, 15._cm, 20._cm, 12._cm));
  SingleDiamondPortalShell polygShell{gctx, testVol};
  // before apply to volueme called - the volume should have zero portals
  BOOST_CHECK_EQUAL(testVol.portals().size(), 0);

  polygShell.applyToVolume();

  // let's see now
  BOOST_CHECK_EQUAL(testVol.portals().size(), 8);
}

BOOST_AUTO_TEST_SUITE_END();

}  // namespace ActsTests
