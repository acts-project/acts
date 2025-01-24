// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/monomorphic/fwd.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/context.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CuboidPortalShell.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <initializer_list>

using namespace Acts::UnitLiterals;

namespace Acts::Test {
GeometryContext gctx;

std::size_t getVolumeIndex() {
  static std::size_t i = 1;
  return i++;
}

auto makeVolume(auto&&... pars) {
  TrackingVolume vol(Transform3::Identity(),
                     std::make_shared<CuboidVolumeBounds>(
                         std::forward<decltype(pars)>(pars)...));
  vol.setVolumeName("cube" + std::to_string(getVolumeIndex()));
  return vol;
};

auto logger = Acts::getDefaultLogger("UnitTests", Acts::Logging::VERBOSE);

BOOST_AUTO_TEST_SUITE(PortalShellTests)

BOOST_AUTO_TEST_CASE(ConstructionFromVolume) {
  auto cube = makeVolume(30_mm, 40_mm, 50_mm);

  TrackingVolume cylVolume(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(10_mm, 20_mm, 10_mm));

  BOOST_CHECK_THROW(SingleCuboidPortalShell{cylVolume}, std::invalid_argument);

  SingleCuboidPortalShell shell1{cube};
  BOOST_CHECK_EQUAL(shell1.size(), 6);

  using enum CuboidVolumeBounds::Face;

  // XY plane
  const auto* pXY = shell1.portal(positiveXYPlane);
  BOOST_REQUIRE_NE(pXY, nullptr);
  BOOST_CHECK_EQUAL(
      pXY->resolveVolume(gctx, Vector3{35_mm, 20_mm, 50_mm}, -Vector3::UnitZ())
          .value(),
      &cube);
  BOOST_CHECK_EQUAL(
      pXY->resolveVolume(gctx, Vector3{35_mm, 20_mm, 50_mm}, Vector3::UnitZ())
          .value(),
      nullptr);

  const auto* nXY = shell1.portal(negativeXYPlane);
  BOOST_REQUIRE_NE(nXY, nullptr);
  BOOST_CHECK_EQUAL(
      nXY->resolveVolume(gctx, Vector3{35_mm, 20_mm, -50_mm}, -Vector3::UnitZ())
          .value(),
      nullptr);
  BOOST_CHECK_EQUAL(
      nXY->resolveVolume(gctx, Vector3{35_mm, 20_mm, -50_mm}, Vector3::UnitZ())
          .value(),
      &cube);

  // YZ plane
  const auto* pYZ = shell1.portal(positiveYZPlane);
  BOOST_REQUIRE_NE(pYZ, nullptr);
  BOOST_CHECK_EQUAL(
      pYZ->resolveVolume(gctx, Vector3{30_mm, 10_mm, 30_mm}, -Vector3::UnitX())
          .value(),
      &cube);
  BOOST_CHECK_EQUAL(
      pYZ->resolveVolume(gctx, Vector3{30_mm, 10_mm, 30_mm}, Vector3::UnitX())
          .value(),
      nullptr);

  const auto* nYZ = shell1.portal(negativeYZPlane);
  BOOST_REQUIRE_NE(nYZ, nullptr);
  BOOST_CHECK_EQUAL(
      nYZ->resolveVolume(gctx, Vector3{-30_mm, 10_mm, 30_mm}, -Vector3::UnitX())
          .value(),
      nullptr);
  BOOST_CHECK_EQUAL(
      nYZ->resolveVolume(gctx, Vector3{-30_mm, 10_mm, 30_mm}, Vector3::UnitX())
          .value(),
      &cube);

  // ZX plane
  const auto* pZX = shell1.portal(positiveZXPlane);
  BOOST_REQUIRE_NE(pZX, nullptr);
  BOOST_CHECK_EQUAL(
      pZX->resolveVolume(gctx, Vector3{15_mm, 40_mm, -10_mm}, -Vector3::UnitY())
          .value(),
      &cube);
  BOOST_CHECK_EQUAL(
      pZX->resolveVolume(gctx, Vector3{15_mm, 40_mm, -10_mm}, Vector3::UnitY())
          .value(),
      nullptr);

  const auto* nZX = shell1.portal(negativeZXPlane);
  BOOST_REQUIRE_NE(nZX, nullptr);
  BOOST_CHECK_EQUAL(
      nZX->resolveVolume(gctx, Vector3{15_mm, 40_mm, -10_mm}, -Vector3::UnitY())
          .value(),
      nullptr);
  BOOST_CHECK_EQUAL(
      nZX->resolveVolume(gctx, Vector3{15_mm, 40_mm, -10_mm}, Vector3::UnitY())
          .value(),
      &cube);
}

BOOST_AUTO_TEST_CASE(PortalAssignment) {
  using enum CuboidVolumeBounds::Face;
  TrackingVolume vol(Transform3::Identity(),
                     std::make_shared<CuboidVolumeBounds>(30_mm, 40_mm, 50_mm));

  SingleCuboidPortalShell shell{vol};

  const auto* pXY = shell.portal(positiveXYPlane);
  const auto* nXY = shell.portal(negativeXYPlane);
  const auto* nYZ = shell.portal(negativeYZPlane);
  const auto* pZX = shell.portal(positiveZXPlane);
  auto* pYZ = shell.portal(positiveYZPlane);
  auto* nZX = shell.portal(negativeZXPlane);

  // Setting new pYZ
  BOOST_REQUIRE_NE(pYZ, nullptr);
  auto* pYZLink = dynamic_cast<const TrivialPortalLink*>(
      pYZ->getLink(Direction::OppositeNormal()));
  BOOST_REQUIRE_NE(pYZLink, nullptr);

  auto grid = pYZLink->makeGrid(AxisDirection::AxisX);

  auto portal2 =
      std::make_shared<Portal>(Direction::OppositeNormal(), std::move(grid));
  shell.setPortal(portal2, positiveYZPlane);
  BOOST_CHECK_EQUAL(shell.portal(positiveYZPlane), portal2.get());

  // Other portals should stay the same
  BOOST_CHECK_EQUAL(shell.portal(positiveXYPlane), pXY);
  BOOST_CHECK_EQUAL(shell.portal(negativeXYPlane), nXY);
  BOOST_CHECK_EQUAL(shell.portal(negativeYZPlane), nYZ);
  BOOST_CHECK_EQUAL(shell.portal(positiveZXPlane), pZX);
  BOOST_CHECK_EQUAL(shell.portal(negativeZXPlane), nZX);

  // Setting new nZX
  BOOST_REQUIRE_NE(nZX, nullptr);
  auto* nZXLink = dynamic_cast<const TrivialPortalLink*>(
      nZX->getLink(Direction::AlongNormal()));
  BOOST_REQUIRE_NE(nZXLink, nullptr);

  grid = nZXLink->makeGrid(AxisDirection::AxisY);

  auto portal3 =
      std::make_shared<Portal>(Direction::AlongNormal(), std::move(grid));
  shell.setPortal(portal3, negativeZXPlane);
  BOOST_CHECK_EQUAL(shell.portal(negativeZXPlane), portal3.get());

  // Other portals should stay the same
  BOOST_CHECK_EQUAL(shell.portal(positiveXYPlane), pXY);
  BOOST_CHECK_EQUAL(shell.portal(negativeXYPlane), nXY);
  BOOST_CHECK_EQUAL(shell.portal(negativeYZPlane), nYZ);
  BOOST_CHECK_EQUAL(shell.portal(positiveZXPlane), pZX);
  BOOST_CHECK_EQUAL(shell.portal(positiveYZPlane), portal2.get());
}

BOOST_AUTO_TEST_SUITE(CuboidStack)
BOOST_DATA_TEST_CASE(XYZDirectionLocal,
                     boost::unit_test::data::make(Acts::AxisDirection::AxisX,
                                                  Acts::AxisDirection::AxisY,
                                                  Acts::AxisDirection::AxisZ),
                     dir) {
  AxisDirection dirOrth1, dirOrth2;
  std::size_t dirIdx;
  switch (dir) {
    case Acts::AxisDirection::AxisX:
      dirOrth1 = Acts::AxisDirection::AxisY;
      dirOrth2 = Acts::AxisDirection::AxisZ;
      dirIdx = 0;
      break;
    case Acts::AxisDirection::AxisY:
      dirOrth1 = Acts::AxisDirection::AxisX;
      dirOrth2 = Acts::AxisDirection::AxisZ;
      dirIdx = 1;
      break;
    case Acts::AxisDirection::AxisZ:
      dirOrth1 = Acts::AxisDirection::AxisX;
      dirOrth2 = Acts::AxisDirection::AxisY;
      dirIdx = 2;
      break;
    default:
      throw std::invalid_argument("Invalid direction");
  }

  auto boundDir = CuboidVolumeBounds::boundsFromAxisDirection(dir);
  auto boundDirOrth1 = CuboidVolumeBounds::boundsFromAxisDirection(dirOrth1);
  auto boundDirOrth2 = CuboidVolumeBounds::boundsFromAxisDirection(dirOrth2);

  auto [frontFace, backFace, sideFaces] =
      CuboidVolumeBounds::facesFromAxisDirection(dir);

  using enum CuboidVolumeBounds::Face;
  auto bounds1 = std::make_shared<CuboidVolumeBounds>(
      std::initializer_list<std::pair<CuboidVolumeBounds::BoundValues, double>>{
          {{boundDir, 100_mm},
           {boundDirOrth1, 30_mm},
           {boundDirOrth2, 100_mm}}});
  auto bounds2 = std::make_shared<CuboidVolumeBounds>(
      std::initializer_list<std::pair<CuboidVolumeBounds::BoundValues, double>>{
          {{boundDir, 100_mm},
           {boundDirOrth1, 30_mm},
           {boundDirOrth2, 100_mm}}});

  TrackingVolume vol1(Transform3{Translation3{Vector3::Unit(dirIdx) * -100_mm}},
                      bounds1);
  TrackingVolume vol2(Transform3{Translation3{Vector3::Unit(dirIdx) * 100_mm}},
                      bounds2);

  SingleCuboidPortalShell shell1{vol1};
  SingleCuboidPortalShell shell2{vol2};

  BOOST_CHECK_NE(shell1.portal(backFace), shell2.portal(frontFace));

  CuboidStackPortalShell stack(gctx, {&shell1, &shell2}, dir, *logger);
  BOOST_CHECK_EQUAL(stack.size(), 6);

  BOOST_CHECK_EQUAL(shell1.portal(frontFace), stack.portal(frontFace));
  BOOST_CHECK_EQUAL(shell1.portal(backFace), shell2.portal(frontFace));
  BOOST_CHECK_EQUAL(shell2.portal(backFace), stack.portal(backFace));

  for (const auto& face : sideFaces) {
    BOOST_CHECK_EQUAL(shell1.portal(face), stack.portal(face));
    BOOST_CHECK_EQUAL(shell2.portal(face), stack.portal(face));
  }

  shell1 = SingleCuboidPortalShell{vol1};
  shell2 = SingleCuboidPortalShell{vol2};

  BOOST_CHECK_THROW(
      CuboidStackPortalShell(gctx, {&shell1, &shell2}, AxisDirection::AxisR),
      std::invalid_argument);
}

BOOST_DATA_TEST_CASE(XYZDirectionGlobal,
                     boost::unit_test::data::xrange(-135, 180, 45), angle) {
  AngleAxis3 baseRotation = AngleAxis3(angle, Vector3::UnitX());
  Vector3 orientation = baseRotation * Vector3::UnitZ();

  auto dir = Acts::AxisDirection::AxisZ;
  auto dirOrth1 = Acts::AxisDirection::AxisX;
  auto dirOrth2 = Acts::AxisDirection::AxisY;
  std::size_t dirIdx = 2;

  auto boundDir = CuboidVolumeBounds::boundsFromAxisDirection(dir);
  auto boundDirOrth1 = CuboidVolumeBounds::boundsFromAxisDirection(dirOrth1);
  auto boundDirOrth2 = CuboidVolumeBounds::boundsFromAxisDirection(dirOrth2);

  auto [frontFace, backFace, sideFaces] =
      CuboidVolumeBounds::facesFromAxisDirection(dir);

  using enum CuboidVolumeBounds::Face;
  auto bounds1 = std::make_shared<CuboidVolumeBounds>(
      std::initializer_list<std::pair<CuboidVolumeBounds::BoundValues, double>>{
          {{boundDir, 100_mm},
           {boundDirOrth1, 30_mm},
           {boundDirOrth2, 100_mm}}});
  auto bounds2 = std::make_shared<CuboidVolumeBounds>(
      std::initializer_list<std::pair<CuboidVolumeBounds::BoundValues, double>>{
          {{boundDir, 100_mm},
           {boundDirOrth1, 30_mm},
           {boundDirOrth2, 100_mm}}});

  Transform3 tranform1 =
      baseRotation * Translation3(Vector3::Unit(dirIdx) * -100_mm);
  TrackingVolume vol1(tranform1, bounds1);
  Transform3 tranform2 =
      baseRotation * Translation3(Vector3::Unit(dirIdx) * 100_mm);
  TrackingVolume vol2(tranform2, bounds2);

  SingleCuboidPortalShell shell1{vol1};
  SingleCuboidPortalShell shell2{vol2};

  BOOST_CHECK_NE(shell1.portal(backFace), shell2.portal(frontFace));

  CuboidStackPortalShell stack(gctx, {&shell1, &shell2}, orientation, *logger);
  BOOST_CHECK_EQUAL(stack.size(), 6);

  BOOST_CHECK_EQUAL(shell1.portal(frontFace), stack.portal(frontFace));
  BOOST_CHECK_EQUAL(shell1.portal(backFace), shell2.portal(frontFace));
  BOOST_CHECK_EQUAL(shell2.portal(backFace), stack.portal(backFace));

  for (const auto& face : sideFaces) {
    BOOST_CHECK_EQUAL(shell1.portal(face), stack.portal(face));
    BOOST_CHECK_EQUAL(shell2.portal(face), stack.portal(face));
  }

  shell1 = SingleCuboidPortalShell{vol1};
  shell2 = SingleCuboidPortalShell{vol2};

  BOOST_CHECK_THROW(
      CuboidStackPortalShell(gctx, {&shell1, &shell2}, AxisDirection::AxisR),
      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(NestedStacks) {
  //   ^
  // z |    +---------------------------------+---------+
  //   |    |                                 |         |
  //   |    |                                 |         |
  //   |    |              vol3               |         |
  //   |    |                                 |         |
  //   |    |                                 |         |
  //   |    +---------------------------------+         |
  //   |    |                                 |         |
  //   |    |                                 |         |
  //   |    |              vol2               |  vol4   |
  //   |    |                                 |         |
  //   |    |                                 |         |
  //   |    +---------------------------------+         |
  //   |    |                                 |         |
  //   |    |                                 |         |
  //   |    |              vol1               |         |
  //   |    |                                 |         |
  //   |    |                                 |         |
  //   |    +---------------------------------+---------+
  //   |
  //   +-------------------------------------------------->
  //                                                      x

  Transform3 base = Transform3::Identity();

  TrackingVolume vol1(
      base, std::make_shared<CuboidVolumeBounds>(30_mm, 100_mm, 200_mm),
      "vol1");

  TrackingVolume vol2(
      base * Translation3(Vector3::UnitZ() * 300_mm),
      std::make_shared<CuboidVolumeBounds>(30_mm, 100_mm, 100_mm), "vol2");

  TrackingVolume vol3(
      base * Translation3(Vector3::UnitZ() * 600_mm),
      std::make_shared<CuboidVolumeBounds>(30_mm, 100_mm, 200_mm), "vol3");

  TrackingVolume vol4(
      base * Translation3{Vector3::UnitX() * 60_mm + Vector3::UnitZ() * 300_mm},
      std::make_shared<CuboidVolumeBounds>(30_mm, 100_mm, 500_mm), "vol4");

  SingleCuboidPortalShell shell1{vol1};
  BOOST_CHECK(shell1.isValid());
  SingleCuboidPortalShell shell2{vol2};
  BOOST_CHECK(shell2.isValid());
  SingleCuboidPortalShell shell3{vol3};
  BOOST_CHECK(shell3.isValid());

  CuboidStackPortalShell stack{
      gctx, {&shell1, &shell2, &shell3}, AxisDirection::AxisZ};

  BOOST_CHECK(stack.isValid());

  SingleCuboidPortalShell shell4{vol4};
  BOOST_CHECK(shell4.isValid());

  CuboidStackPortalShell stack2{
      gctx, {&stack, &shell4}, AxisDirection::AxisX, *logger};
  BOOST_CHECK(stack2.isValid());

  using enum CuboidVolumeBounds::Face;

  auto lookup = [](auto& shell, CuboidPortalShell::Face face, Vector3 position,
                   Vector3 direction) -> const TrackingVolume* {
    const auto* portal = shell.portal(face);
    BOOST_REQUIRE_NE(portal, nullptr);
    const auto* vol = portal->resolveVolume(gctx, position, direction).value();
    if (vol != nullptr) {
      std::cout << "PORTAL " << portal->surface().center(gctx).transpose()
                << " -- VOLUME: " << vol->volumeName() << "\n";
    } else {
      std::cout << "PORTAL " << portal->surface().center(gctx).transpose()
                << " -- VOLUME: NULL\n";
    }
    return vol;
  };

  // Shell 1
  BOOST_CHECK_EQUAL(lookup(shell1, negativeXYPlane,
                           Vector3(10_mm, 20_mm, -200_mm), -Vector3::UnitZ()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(shell1, negativeXYPlane,
                           Vector3(10_mm, 20_mm, -200_mm), Vector3::UnitZ()),
                    &vol1);

  BOOST_CHECK_EQUAL(lookup(shell1, positiveXYPlane,
                           Vector3(10_mm, 20_mm, 200_mm), -Vector3::UnitZ()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(shell1, positiveXYPlane,
                           Vector3(10_mm, 20_mm, 200_mm), Vector3::UnitZ()),
                    &vol2);

  BOOST_CHECK_EQUAL(lookup(shell1, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 20_mm), -Vector3::UnitX()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(shell1, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 20_mm), Vector3::UnitX()),
                    &vol1);

  BOOST_CHECK_EQUAL(lookup(shell1, positiveYZPlane,
                           Vector3(30_mm, 10_mm, 20_mm), -Vector3::UnitX()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(shell1, positiveYZPlane,
                           Vector3(30_mm, 10_mm, 20_mm), Vector3::UnitX()),
                    &vol4);

  BOOST_CHECK_EQUAL(lookup(shell1, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 20_mm), -Vector3::UnitY()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(shell1, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 20_mm), Vector3::UnitY()),
                    &vol1);

  BOOST_CHECK_EQUAL(lookup(shell1, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 20_mm), -Vector3::UnitY()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(shell1, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 20_mm), Vector3::UnitY()),
                    nullptr);

  // Volume 2
  BOOST_CHECK_EQUAL(lookup(shell2, negativeXYPlane,
                           Vector3(10_mm, 20_mm, 200_mm), -Vector3::UnitZ()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(shell2, negativeXYPlane,
                           Vector3(10_mm, 20_mm, 200_mm), Vector3::UnitZ()),
                    &vol2);

  BOOST_CHECK_EQUAL(lookup(shell2, positiveXYPlane,
                           Vector3(10_mm, 20_mm, 400_mm), -Vector3::UnitZ()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(shell2, positiveXYPlane,
                           Vector3(10_mm, 20_mm, 400_mm), Vector3::UnitZ()),
                    &vol3);

  BOOST_CHECK_EQUAL(lookup(shell2, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 220_mm), -Vector3::UnitX()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(shell2, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 220_mm), Vector3::UnitX()),
                    &vol2);

  BOOST_CHECK_EQUAL(lookup(shell2, positiveYZPlane,
                           Vector3(30_mm, 10_mm, 220_mm), -Vector3::UnitX()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(shell2, positiveYZPlane,
                           Vector3(30_mm, 10_mm, 220_mm), Vector3::UnitX()),
                    &vol4);

  BOOST_CHECK_EQUAL(lookup(shell2, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 220_mm), -Vector3::UnitY()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(shell2, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 220_mm), Vector3::UnitY()),
                    &vol2);

  BOOST_CHECK_EQUAL(lookup(shell2, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 220_mm), -Vector3::UnitY()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(shell2, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 220_mm), Vector3::UnitY()),
                    nullptr);

  // Volume 3
  BOOST_CHECK_EQUAL(lookup(shell3, negativeXYPlane,
                           Vector3(10_mm, 20_mm, 400_mm), -Vector3::UnitZ()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(shell3, negativeXYPlane,
                           Vector3(10_mm, 20_mm, 400_mm), Vector3::UnitZ()),
                    &vol3);

  BOOST_CHECK_EQUAL(lookup(shell3, positiveXYPlane,
                           Vector3(10_mm, 20_mm, 800_mm), -Vector3::UnitZ()),
                    &vol3);
  BOOST_CHECK_EQUAL(lookup(shell3, positiveXYPlane,
                           Vector3(10_mm, 20_mm, 800_mm), Vector3::UnitZ()),
                    nullptr);

  BOOST_CHECK_EQUAL(lookup(shell3, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 420_mm), -Vector3::UnitX()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(shell3, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 420_mm), Vector3::UnitX()),
                    &vol3);

  BOOST_CHECK_EQUAL(lookup(shell3, positiveYZPlane,
                           Vector3(30_mm, 10_mm, 420_mm), -Vector3::UnitX()),
                    &vol3);
  BOOST_CHECK_EQUAL(lookup(shell3, positiveYZPlane,
                           Vector3(30_mm, 10_mm, 420_mm), Vector3::UnitX()),
                    &vol4);

  BOOST_CHECK_EQUAL(lookup(shell3, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 420_mm), -Vector3::UnitY()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(shell3, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 420_mm), Vector3::UnitY()),
                    &vol3);

  BOOST_CHECK_EQUAL(lookup(shell3, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 420_mm), -Vector3::UnitY()),
                    &vol3);
  BOOST_CHECK_EQUAL(lookup(shell3, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 420_mm), Vector3::UnitY()),
                    nullptr);

  // Volume 4

  for (CuboidVolumeBounds::Face face :
       {negativeXYPlane, positiveXYPlane, negativeYZPlane, positiveYZPlane,
        negativeZXPlane, positiveZXPlane}) {
    const auto& portalAtFace = shell4.portal(face);
    if (portalAtFace != nullptr) {
      std::cout << "FACE: " << face << " -- "
                << portalAtFace->surface().center(gctx).transpose() << "\n";
    } else {
      std::cout << "FACE: " << face << " NULL\n";
    }
  }

  BOOST_CHECK_EQUAL(lookup(shell4, negativeXYPlane,
                           Vector3(50_mm, 20_mm, -200_mm), -Vector3::UnitZ()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(shell4, negativeXYPlane,
                           Vector3(50_mm, 20_mm, -200_mm), Vector3::UnitZ()),
                    &vol4);

  BOOST_CHECK_EQUAL(lookup(shell4, positiveXYPlane,
                           Vector3(50_mm, 20_mm, 800_mm), -Vector3::UnitZ()),
                    &vol4);
  BOOST_CHECK_EQUAL(lookup(shell4, positiveXYPlane,
                           Vector3(50_mm, 20_mm, 800_mm), Vector3::UnitZ()),
                    nullptr);

  BOOST_CHECK_EQUAL(lookup(shell4, negativeYZPlane, Vector3(30_mm, 10_mm, 0_mm),
                           -Vector3::UnitX()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(shell4, negativeYZPlane,
                           Vector3(30_mm, 10_mm, 220_mm), -Vector3::UnitX()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(shell4, negativeYZPlane,
                           Vector3(30_mm, 10_mm, 420_mm), -Vector3::UnitX()),
                    &vol3);
  BOOST_CHECK_EQUAL(lookup(shell4, negativeYZPlane,
                           Vector3(30_mm, 10_mm, 300_mm), Vector3::UnitX()),
                    &vol4);

  BOOST_CHECK_EQUAL(lookup(shell4, positiveYZPlane,
                           Vector3(60_mm, 10_mm, 300_mm), -Vector3::UnitX()),
                    &vol4);
  BOOST_CHECK_EQUAL(lookup(shell4, positiveYZPlane,
                           Vector3(60_mm, 10_mm, 300_mm), Vector3::UnitX()),
                    nullptr);

  BOOST_CHECK_EQUAL(lookup(shell4, negativeZXPlane,
                           Vector3(50_mm, -100_mm, 300_mm), -Vector3::UnitY()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(shell4, negativeZXPlane,
                           Vector3(50_mm, -100_mm, 300_mm), Vector3::UnitY()),
                    &vol4);

  BOOST_CHECK_EQUAL(lookup(shell4, positiveZXPlane,
                           Vector3(50_mm, 100_mm, 300_mm), -Vector3::UnitY()),
                    &vol4);
  BOOST_CHECK_EQUAL(lookup(shell4, positiveZXPlane,
                           Vector3(50_mm, 100_mm, 300_mm), Vector3::UnitY()),
                    nullptr);

  // Stack
  BOOST_CHECK_EQUAL(lookup(stack, negativeXYPlane,
                           Vector3(10_mm, 20_mm, -200_mm), -Vector3::UnitZ()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(stack, negativeXYPlane,
                           Vector3(10_mm, 20_mm, -200_mm), Vector3::UnitZ()),
                    &vol1);

  BOOST_CHECK_EQUAL(lookup(stack, positiveXYPlane,
                           Vector3(10_mm, 20_mm, 800_mm), -Vector3::UnitZ()),
                    &vol3);
  BOOST_CHECK_EQUAL(lookup(stack, positiveXYPlane,
                           Vector3(10_mm, 20_mm, 800_mm), Vector3::UnitZ()),
                    nullptr);

  BOOST_CHECK_EQUAL(lookup(stack, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 300_mm), -Vector3::UnitX()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(stack, negativeYZPlane, Vector3(-30_mm, 10_mm, 0_mm),
                           Vector3::UnitX()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(stack, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 220_mm), Vector3::UnitX()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(stack, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 420_mm), Vector3::UnitX()),
                    &vol3);

  BOOST_CHECK_EQUAL(lookup(stack, positiveYZPlane, Vector3(30_mm, 10_mm, 0_mm),
                           -Vector3::UnitX()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(stack, positiveYZPlane,
                           Vector3(30_mm, 10_mm, 220_mm), -Vector3::UnitX()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(stack, positiveYZPlane,
                           Vector3(30_mm, 10_mm, 420_mm), -Vector3::UnitX()),
                    &vol3);
  BOOST_CHECK_EQUAL(lookup(stack, positiveYZPlane,
                           Vector3(30_mm, 10_mm, 300_mm), Vector3::UnitX()),
                    &vol4);

  BOOST_CHECK_EQUAL(lookup(stack, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 300_mm), -Vector3::UnitY()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(stack, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 0_mm), Vector3::UnitY()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(stack, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 220_mm), Vector3::UnitY()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(stack, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 420_mm), Vector3::UnitY()),
                    &vol3);

  BOOST_CHECK_EQUAL(lookup(stack, positiveZXPlane, Vector3(10_mm, 100_mm, 0_mm),
                           -Vector3::UnitY()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(stack, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 220_mm), -Vector3::UnitY()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(stack, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 420_mm), -Vector3::UnitY()),
                    &vol3);
  BOOST_CHECK_EQUAL(lookup(stack, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 300_mm), Vector3::UnitY()),
                    nullptr);

  // Stack 2
  BOOST_CHECK_EQUAL(lookup(stack2, negativeXYPlane,
                           Vector3(10_mm, 20_mm, -200_mm), -Vector3::UnitZ()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(stack2, negativeXYPlane,
                           Vector3(10_mm, 20_mm, -200_mm), Vector3::UnitZ()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(stack2, negativeXYPlane,
                           Vector3(50_mm, 20_mm, -200_mm), Vector3::UnitZ()),
                    &vol4);

  BOOST_CHECK_EQUAL(lookup(stack2, positiveXYPlane,
                           Vector3(10_mm, 20_mm, 800_mm), -Vector3::UnitZ()),
                    &vol3);
  BOOST_CHECK_EQUAL(lookup(stack2, positiveXYPlane,
                           Vector3(50_mm, 20_mm, 800_mm), -Vector3::UnitZ()),
                    &vol4);
  BOOST_CHECK_EQUAL(lookup(stack2, positiveXYPlane,
                           Vector3(10_mm, 20_mm, 800_mm), Vector3::UnitZ()),
                    nullptr);

  BOOST_CHECK_EQUAL(lookup(stack2, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 300_mm), -Vector3::UnitX()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(stack2, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 0_mm), Vector3::UnitX()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(stack2, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 220_mm), Vector3::UnitX()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(stack2, negativeYZPlane,
                           Vector3(-30_mm, 10_mm, 420_mm), Vector3::UnitX()),
                    &vol3);

  BOOST_CHECK_EQUAL(lookup(shell4, positiveYZPlane,
                           Vector3(60_mm, 10_mm, 300_mm), -Vector3::UnitX()),
                    &vol4);
  BOOST_CHECK_EQUAL(lookup(shell4, positiveYZPlane,
                           Vector3(60_mm, 10_mm, 300_mm), Vector3::UnitX()),
                    nullptr);

  BOOST_CHECK_EQUAL(lookup(stack2, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 300_mm), -Vector3::UnitY()),
                    nullptr);
  BOOST_CHECK_EQUAL(lookup(stack2, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 0_mm), Vector3::UnitY()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(stack2, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 220_mm), Vector3::UnitY()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(stack2, negativeZXPlane,
                           Vector3(10_mm, -100_mm, 420_mm), Vector3::UnitY()),
                    &vol3);
  BOOST_CHECK_EQUAL(lookup(stack2, negativeZXPlane,
                           Vector3(50_mm, -100_mm, 420_mm), Vector3::UnitY()),
                    &vol4);

  BOOST_CHECK_EQUAL(lookup(stack2, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 0_mm), -Vector3::UnitY()),
                    &vol1);
  BOOST_CHECK_EQUAL(lookup(stack2, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 220_mm), -Vector3::UnitY()),
                    &vol2);
  BOOST_CHECK_EQUAL(lookup(stack2, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 420_mm), -Vector3::UnitY()),
                    &vol3);
  BOOST_CHECK_EQUAL(lookup(stack2, positiveZXPlane,
                           Vector3(50_mm, 100_mm, 420_mm), -Vector3::UnitY()),
                    &vol4);
  BOOST_CHECK_EQUAL(lookup(stack2, positiveZXPlane,
                           Vector3(10_mm, 100_mm, 300_mm), Vector3::UnitY()),
                    nullptr);
}

BOOST_AUTO_TEST_CASE(Fill) {
  Transform3 base = Transform3::Identity();

  TrackingVolume vol1(
      base, std::make_shared<CuboidVolumeBounds>(30_mm, 100_mm, 200_mm),
      "vol1");

  TrackingVolume vol2(
      base * Translation3(Vector3::UnitZ() * 300_mm),
      std::make_shared<CuboidVolumeBounds>(30_mm, 100_mm, 100_mm), "vol2");

  SingleCuboidPortalShell shell{vol1};

  using enum CuboidVolumeBounds::Face;

  BOOST_CHECK_EQUAL(
      shell.portal(positiveXYPlane)->getLink(Direction::AlongNormal()),
      nullptr);

  shell.fill(vol2);

  BOOST_CHECK_NE(
      shell.portal(positiveXYPlane)->getLink(Direction::AlongNormal()),
      nullptr);
}

BOOST_AUTO_TEST_CASE(RegisterInto) {
  using enum CuboidVolumeBounds::Face;
  TrackingVolume vol1(
      Transform3::Identity(),
      std::make_shared<CuboidVolumeBounds>(100_mm, 100_mm, 100_mm));

  SingleCuboidPortalShell shell{vol1};

  BOOST_CHECK_EQUAL(vol1.portals().size(), 0);

  shell.applyToVolume();
  BOOST_CHECK_EQUAL(vol1.portals().size(), 6);
}

BOOST_AUTO_TEST_SUITE_END()  // CuboidStack
BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
