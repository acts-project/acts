// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/context.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalShell.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Surfaces/SurfaceMergingException.hpp"

using namespace Acts::UnitLiterals;

namespace Acts::Test {
GeometryContext gctx;

BOOST_AUTO_TEST_SUITE(PortalShellTests)

BOOST_AUTO_TEST_CASE(ConstructionFromVolume) {
  // - Cylinder
  // |           | no phi | phi |
  // | --------- | ------ | --- |
  // | rMin > 0  | 1      | 3   |
  // | rMin == 0 | 2      | 4   |

  std::size_t i = 1;
  auto makeVolume = [&](auto&&... pars) {
    TrackingVolume vol(Transform3::Identity(),
                       std::make_shared<CylinderVolumeBounds>(
                           std::forward<decltype(pars)>(pars)...));
    vol.setVolumeName("cyl" + std::to_string(i++));
    return std::move(vol);
  };

  auto cyl1 = makeVolume(30_mm, 40_mm, 100_mm);
  auto cyl2 = makeVolume(0_mm, 40_mm, 100_mm);
  auto cyl3 = makeVolume(30_mm, 40_mm, 100_mm, 45_degree);
  auto cyl4 = makeVolume(0_mm, 40_mm, 100_mm, 45_degree);

  SingleCylinderPortalShell shell1{gctx, cyl1};
  BOOST_CHECK_EQUAL(shell1.size(), 4);

  using enum CylinderPortalShell::Face;

  const auto* pDisc = shell1.portal(PositiveDisc);
  BOOST_REQUIRE_NE(pDisc, nullptr);
  BOOST_CHECK_EQUAL(pDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, 100_mm},
                                         -Vector3::UnitZ()),
                    &cyl1);
  BOOST_CHECK_EQUAL(pDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, 100_mm},
                                         Vector3::UnitZ()),
                    nullptr);

  const auto* nDisc = shell1.portal(NegativeDisc);
  BOOST_REQUIRE_NE(nDisc, nullptr);
  BOOST_CHECK_EQUAL(nDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, -100_mm},
                                         -Vector3::UnitZ()),
                    nullptr);
  BOOST_CHECK_EQUAL(nDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, -100_mm},
                                         Vector3::UnitZ()),
                    &cyl1);

  const auto* oCyl = shell1.portal(OuterCylinder);
  BOOST_REQUIRE_NE(oCyl, nullptr);
  BOOST_CHECK_EQUAL(
      oCyl->resolveVolume(gctx, Vector3{40_mm, 0_mm, 10_mm}, Vector3::UnitX()),
      nullptr);
  BOOST_CHECK_EQUAL(
      oCyl->resolveVolume(gctx, Vector3{40_mm, 0_mm, 10_mm}, -Vector3::UnitX()),
      &cyl1);

  const auto* iCyl = shell1.portal(InnerCylinder);
  BOOST_REQUIRE_NE(iCyl, nullptr);
  BOOST_CHECK_EQUAL(
      iCyl->resolveVolume(gctx, Vector3{30_mm, 0_mm, 10_mm}, Vector3::UnitX()),
      &cyl1);
  BOOST_CHECK_EQUAL(
      iCyl->resolveVolume(gctx, Vector3{30_mm, 0_mm, 10_mm}, -Vector3::UnitX()),
      nullptr);

  SingleCylinderPortalShell shell2{gctx, cyl2};
  BOOST_CHECK_EQUAL(shell2.size(), 3);

  pDisc = shell2.portal(PositiveDisc);
  BOOST_REQUIRE_NE(pDisc, nullptr);
  BOOST_CHECK_EQUAL(pDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, 100_mm},
                                         -Vector3::UnitZ()),
                    &cyl2);
  BOOST_CHECK_EQUAL(pDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, 100_mm},
                                         Vector3::UnitZ()),
                    nullptr);

  nDisc = shell2.portal(NegativeDisc);
  BOOST_REQUIRE_NE(nDisc, nullptr);
  BOOST_CHECK_EQUAL(nDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, -100_mm},
                                         -Vector3::UnitZ()),
                    nullptr);
  BOOST_CHECK_EQUAL(nDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, -100_mm},
                                         Vector3::UnitZ()),
                    &cyl2);

  oCyl = shell2.portal(OuterCylinder);
  BOOST_REQUIRE_NE(oCyl, nullptr);
  BOOST_CHECK_EQUAL(
      oCyl->resolveVolume(gctx, Vector3{40_mm, 0_mm, 10_mm}, Vector3::UnitX()),
      nullptr);
  BOOST_CHECK_EQUAL(
      oCyl->resolveVolume(gctx, Vector3{40_mm, 0_mm, 10_mm}, -Vector3::UnitX()),
      &cyl2);

  iCyl = shell2.portal(InnerCylinder);
  BOOST_CHECK_EQUAL(iCyl, nullptr);

  SingleCylinderPortalShell shell3{gctx, cyl3};
  BOOST_CHECK_EQUAL(shell3.size(), 6);

  pDisc = shell3.portal(PositiveDisc);
  BOOST_REQUIRE_NE(pDisc, nullptr);
  BOOST_CHECK_EQUAL(pDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, 100_mm},
                                         -Vector3::UnitZ()),
                    &cyl3);
  BOOST_CHECK_EQUAL(pDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, 100_mm},
                                         Vector3::UnitZ()),
                    nullptr);

  nDisc = shell3.portal(NegativeDisc);
  BOOST_REQUIRE_NE(nDisc, nullptr);
  BOOST_CHECK_EQUAL(nDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, -100_mm},
                                         -Vector3::UnitZ()),
                    nullptr);
  BOOST_CHECK_EQUAL(nDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, -100_mm},
                                         Vector3::UnitZ()),
                    &cyl3);

  oCyl = shell3.portal(OuterCylinder);
  BOOST_REQUIRE_NE(oCyl, nullptr);
  BOOST_CHECK_EQUAL(
      oCyl->resolveVolume(gctx, Vector3{40_mm, 0_mm, 10_mm}, Vector3::UnitX()),
      nullptr);
  BOOST_CHECK_EQUAL(
      oCyl->resolveVolume(gctx, Vector3{40_mm, 0_mm, 10_mm}, -Vector3::UnitX()),
      &cyl3);

  iCyl = shell3.portal(InnerCylinder);
  BOOST_REQUIRE_NE(iCyl, nullptr);
  BOOST_CHECK_EQUAL(
      iCyl->resolveVolume(gctx, Vector3{30_mm, 0_mm, 10_mm}, Vector3::UnitX()),
      &cyl3);
  BOOST_CHECK_EQUAL(
      iCyl->resolveVolume(gctx, Vector3{30_mm, 0_mm, 10_mm}, -Vector3::UnitX()),
      nullptr);

  auto anglePoint = [](double angle, double r, double z) {
    return Vector3{r * std::cos(angle), r * std::sin(angle), z};
  };

  const auto* nPhi = shell3.portal(NegativePhiPlane);
  BOOST_REQUIRE_NE(nPhi, nullptr);
  Vector3 point = anglePoint(-45_degree, 35_mm, 10_mm);
  Vector3 dir = Vector3::UnitZ().cross(point).normalized();
  Vector3 idir = (-Vector3::UnitZ()).cross(point).normalized();
  BOOST_CHECK_EQUAL(nPhi->resolveVolume(gctx, point, dir), nullptr);
  BOOST_CHECK_EQUAL(nPhi->resolveVolume(gctx, point, idir), &cyl3);

  const auto* pPhi = shell3.portal(PositivePhiPlane);
  BOOST_REQUIRE_NE(pPhi, nullptr);
  point = anglePoint(45_degree, 35_mm, 10_mm);
  dir = Vector3::UnitZ().cross(point).normalized();
  idir = (-Vector3::UnitZ()).cross(point).normalized();
  BOOST_CHECK_EQUAL(pPhi->resolveVolume(gctx, point, dir), nullptr);
  BOOST_CHECK_EQUAL(pPhi->resolveVolume(gctx, point, idir), &cyl3);

  SingleCylinderPortalShell shell4{gctx, cyl4};
  BOOST_CHECK_EQUAL(shell4.size(), 5);

  pDisc = shell4.portal(PositiveDisc);
  BOOST_REQUIRE_NE(pDisc, nullptr);
  BOOST_CHECK_EQUAL(pDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, 100_mm},
                                         -Vector3::UnitZ()),
                    &cyl4);
  BOOST_CHECK_EQUAL(pDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, 100_mm},
                                         Vector3::UnitZ()),
                    nullptr);

  nDisc = shell4.portal(NegativeDisc);
  BOOST_REQUIRE_NE(nDisc, nullptr);
  BOOST_CHECK_EQUAL(nDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, -100_mm},
                                         -Vector3::UnitZ()),
                    nullptr);
  BOOST_CHECK_EQUAL(nDisc->resolveVolume(gctx, Vector3{35_mm, 0_mm, -100_mm},
                                         Vector3::UnitZ()),
                    &cyl4);

  oCyl = shell4.portal(OuterCylinder);
  BOOST_REQUIRE_NE(oCyl, nullptr);
  BOOST_CHECK_EQUAL(
      oCyl->resolveVolume(gctx, Vector3{40_mm, 0_mm, 10_mm}, Vector3::UnitX()),
      nullptr);
  BOOST_CHECK_EQUAL(
      oCyl->resolveVolume(gctx, Vector3{40_mm, 0_mm, 10_mm}, -Vector3::UnitX()),
      &cyl4);

  iCyl = shell4.portal(InnerCylinder);
  BOOST_CHECK_EQUAL(iCyl, nullptr);

  nPhi = shell4.portal(NegativePhiPlane);
  BOOST_REQUIRE_NE(nPhi, nullptr);
  point = anglePoint(-45_degree, 35_mm, 10_mm);
  dir = Vector3::UnitZ().cross(point).normalized();
  idir = (-Vector3::UnitZ()).cross(point).normalized();
  BOOST_CHECK_EQUAL(nPhi->resolveVolume(gctx, point, dir), nullptr);
  BOOST_CHECK_EQUAL(nPhi->resolveVolume(gctx, point, idir), &cyl4);

  pPhi = shell4.portal(PositivePhiPlane);
  BOOST_REQUIRE_NE(pPhi, nullptr);
  point = anglePoint(45_degree, 35_mm, 10_mm);
  dir = Vector3::UnitZ().cross(point).normalized();
  idir = (-Vector3::UnitZ()).cross(point).normalized();
  BOOST_CHECK_EQUAL(pPhi->resolveVolume(gctx, point, dir), nullptr);
  BOOST_CHECK_EQUAL(pPhi->resolveVolume(gctx, point, idir), &cyl4);
}

//              outer cylinder
//          +----------+----------+
//          |          |          |
// negative |          v          | positive
//     disc +-->               <--+ disc
//          |          ^          |
//          |          |          |
//          +----------+----------+
//              inner cylinder

BOOST_AUTO_TEST_CASE(PortalAssignment) {
  using enum CylinderPortalShell::Face;
  TrackingVolume vol(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(30_mm, 100_mm, 100_mm));

  SingleCylinderPortalShell shell{gctx, vol};

  const auto* iCyl = shell.portal(InnerCylinder);
  const auto* pDisc = shell.portal(PositiveDisc);
  auto* oCyl = shell.portal(OuterCylinder);
  auto* nDisc = shell.portal(NegativeDisc);

  // Setting new outer cylinder
  BOOST_REQUIRE_NE(oCyl, nullptr);
  auto* oCylLink = dynamic_cast<const TrivialPortalLink*>(
      oCyl->getLink(Direction::OppositeNormal));
  BOOST_REQUIRE_NE(oCylLink, nullptr);

  auto grid = oCylLink->makeGrid(BinningValue::binZ);

  auto portal2 = std::make_shared<Portal>(gctx, Direction::OppositeNormal,
                                          std::move(grid));
  shell.setPortal(portal2, OuterCylinder);
  BOOST_CHECK_EQUAL(shell.portal(OuterCylinder), portal2.get());

  // Other portals should stay the same
  BOOST_CHECK_EQUAL(shell.portal(InnerCylinder), iCyl);
  BOOST_CHECK_EQUAL(shell.portal(PositiveDisc), pDisc);
  BOOST_CHECK_EQUAL(shell.portal(NegativeDisc), nDisc);

  // Setting new negative disc
  BOOST_REQUIRE_NE(nDisc, nullptr);
  auto* nDiscLink = dynamic_cast<const TrivialPortalLink*>(
      nDisc->getLink(Direction::AlongNormal));
  BOOST_REQUIRE_NE(nDiscLink, nullptr);

  grid = nDiscLink->makeGrid(BinningValue::binR);

  auto portal3 =
      std::make_shared<Portal>(gctx, Direction::AlongNormal, std::move(grid));
  shell.setPortal(portal3, NegativeDisc);
  BOOST_CHECK_EQUAL(shell.portal(NegativeDisc), portal3.get());

  // Other portals should stay the same
  BOOST_CHECK_EQUAL(shell.portal(OuterCylinder), portal2.get());
  BOOST_CHECK_EQUAL(shell.portal(InnerCylinder), iCyl);
  BOOST_CHECK_EQUAL(shell.portal(PositiveDisc), pDisc);
}

BOOST_AUTO_TEST_SUITE(CylinderStack)
BOOST_AUTO_TEST_CASE(ZDirection) {
  using enum CylinderPortalShell::Face;
  BOOST_TEST_CONTEXT("rMin>0") {
    TrackingVolume vol1(
        Transform3{Translation3{Vector3::UnitZ() * -100_mm}},
        std::make_shared<CylinderVolumeBounds>(30_mm, 100_mm, 100_mm));

    TrackingVolume vol2(
        Transform3{Translation3{Vector3::UnitZ() * 100_mm}},
        std::make_shared<CylinderVolumeBounds>(30_mm, 100_mm, 100_mm));

    SingleCylinderPortalShell shell1{gctx, vol1};
    SingleCylinderPortalShell shell2{gctx, vol2};

    BOOST_CHECK_NE(shell1.portal(PositiveDisc), shell2.portal(NegativeDisc));

    CylinderStackPortalShell stack{
        gctx, {&shell1, &shell2}, BinningValue::binZ};
    BOOST_CHECK_EQUAL(stack.size(), 4);

    const auto* iCyl = stack.portal(InnerCylinder);
    BOOST_CHECK_EQUAL(shell1.portal(InnerCylinder), iCyl);
    BOOST_CHECK_EQUAL(shell2.portal(InnerCylinder), iCyl);

    const auto* oCyl = stack.portal(OuterCylinder);
    BOOST_CHECK_EQUAL(shell1.portal(OuterCylinder), oCyl);
    BOOST_CHECK_EQUAL(shell2.portal(OuterCylinder), oCyl);

    BOOST_CHECK_EQUAL(stack.portal(PositiveDisc), shell2.portal(PositiveDisc));
    BOOST_CHECK_EQUAL(stack.portal(NegativeDisc), shell1.portal(NegativeDisc));

    // Disc portals have been fused
    BOOST_CHECK_EQUAL(shell1.portal(PositiveDisc), shell2.portal(NegativeDisc));

    shell1 = SingleCylinderPortalShell{gctx, vol1};
    shell2 = SingleCylinderPortalShell{gctx, vol2};

    BOOST_CHECK_THROW(
        CylinderStackPortalShell(gctx, {&shell1}, BinningValue::binZ),
        std::invalid_argument);

    BOOST_CHECK_THROW(
        CylinderStackPortalShell(gctx, {&shell1, &shell2}, BinningValue::binR),
        SurfaceMergingException);
  }

  BOOST_TEST_CONTEXT("rMin==0") {
    TrackingVolume vol1(
        Transform3{Translation3{Vector3::UnitZ() * -100_mm}},
        std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm));

    TrackingVolume vol2(
        Transform3{Translation3{Vector3::UnitZ() * 100_mm}},
        std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm));

    SingleCylinderPortalShell shell1{gctx, vol1};
    SingleCylinderPortalShell shell2{gctx, vol2};

    BOOST_CHECK_EQUAL(shell1.portal(InnerCylinder), nullptr);
    BOOST_CHECK_EQUAL(shell2.portal(InnerCylinder), nullptr);

    BOOST_CHECK_NE(shell1.portal(PositiveDisc), shell2.portal(NegativeDisc));

    CylinderStackPortalShell stack{
        gctx, {&shell1, &shell2}, BinningValue::binZ};
    BOOST_CHECK_EQUAL(stack.size(), 3);

    // Disc portals have been fused
    BOOST_CHECK_EQUAL(shell1.portal(PositiveDisc), shell2.portal(NegativeDisc));

    const auto* iCyl = stack.portal(InnerCylinder);
    BOOST_CHECK_EQUAL(iCyl, nullptr);

    const auto* oCyl = stack.portal(OuterCylinder);
    BOOST_CHECK_EQUAL(shell1.portal(OuterCylinder), oCyl);
    BOOST_CHECK_EQUAL(shell2.portal(OuterCylinder), oCyl);

    BOOST_CHECK_EQUAL(stack.portal(PositiveDisc), shell2.portal(PositiveDisc));
    BOOST_CHECK_EQUAL(stack.portal(NegativeDisc), shell1.portal(NegativeDisc));

    shell1 = SingleCylinderPortalShell{gctx, vol1};
    shell2 = SingleCylinderPortalShell{gctx, vol2};

    BOOST_CHECK_THROW(
        CylinderStackPortalShell(gctx, {&shell1}, BinningValue::binZ),
        std::invalid_argument);

    BOOST_CHECK_THROW(
        CylinderStackPortalShell(gctx, {&shell1, &shell2}, BinningValue::binR),
        SurfaceMergingException);
  }
}

BOOST_AUTO_TEST_CASE(RDirection) {
  using enum CylinderPortalShell::Face;
  BOOST_TEST_CONTEXT("rMin>0") {
    TrackingVolume vol1(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(30_mm, 100_mm, 100_mm));

    TrackingVolume vol2(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(100_mm, 150_mm, 100_mm));

    SingleCylinderPortalShell shell1{gctx, vol1};
    SingleCylinderPortalShell shell2{gctx, vol2};

    BOOST_CHECK_NE(shell1.portal(OuterCylinder), shell2.portal(InnerCylinder));

    CylinderStackPortalShell stack{
        gctx, {&shell1, &shell2}, BinningValue::binR};
    BOOST_CHECK_EQUAL(stack.size(), 4);

    // Internal cylinder portals have been fused
    BOOST_CHECK_EQUAL(shell1.portal(OuterCylinder),
                      shell2.portal(InnerCylinder));

    const auto* nDisc = stack.portal(NegativeDisc);
    const auto* pDisc = stack.portal(PositiveDisc);

    BOOST_CHECK_EQUAL(shell1.portal(NegativeDisc), nDisc);
    BOOST_CHECK_EQUAL(shell2.portal(NegativeDisc), nDisc);
    BOOST_CHECK_EQUAL(shell1.portal(PositiveDisc), pDisc);
    BOOST_CHECK_EQUAL(shell2.portal(PositiveDisc), pDisc);

    BOOST_CHECK_EQUAL(stack.portal(InnerCylinder),
                      shell1.portal(InnerCylinder));
    BOOST_CHECK_EQUAL(stack.portal(OuterCylinder),
                      shell2.portal(OuterCylinder));

    shell1 = SingleCylinderPortalShell{gctx, vol1};
    shell2 = SingleCylinderPortalShell{gctx, vol2};

    BOOST_CHECK_THROW(
        CylinderStackPortalShell(gctx, {&shell1}, BinningValue::binR),
        std::invalid_argument);

    BOOST_CHECK_THROW(
        CylinderStackPortalShell(gctx, {&shell1, &shell2}, BinningValue::binZ),
        SurfaceMergingException);
  }

  BOOST_TEST_CONTEXT("rMin==0") {
    TrackingVolume vol1(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(0_mm, 100_mm, 100_mm));

    TrackingVolume vol2(
        Transform3::Identity(),
        std::make_shared<CylinderVolumeBounds>(100_mm, 150_mm, 100_mm));

    SingleCylinderPortalShell shell1{gctx, vol1};
    SingleCylinderPortalShell shell2{gctx, vol2};

    BOOST_CHECK_EQUAL(shell1.portal(InnerCylinder), nullptr);
    BOOST_CHECK_NE(shell1.portal(OuterCylinder), shell2.portal(InnerCylinder));

    CylinderStackPortalShell stack{
        gctx, {&shell1, &shell2}, BinningValue::binR};
    BOOST_CHECK_EQUAL(stack.size(), 4);

    // Internal cylinder portals have been fused
    BOOST_CHECK_EQUAL(shell1.portal(OuterCylinder),
                      shell2.portal(InnerCylinder));

    const auto* nDisc = stack.portal(NegativeDisc);
    const auto* pDisc = stack.portal(PositiveDisc);

    BOOST_CHECK_EQUAL(shell1.portal(NegativeDisc), nDisc);
    BOOST_CHECK_EQUAL(shell2.portal(NegativeDisc), nDisc);
    BOOST_CHECK_EQUAL(shell1.portal(PositiveDisc), pDisc);
    BOOST_CHECK_EQUAL(shell2.portal(PositiveDisc), pDisc);

    BOOST_CHECK_EQUAL(stack.portal(InnerCylinder), nullptr);
    BOOST_CHECK_EQUAL(stack.portal(OuterCylinder),
                      shell2.portal(OuterCylinder));

    shell1 = SingleCylinderPortalShell{gctx, vol1};
    shell2 = SingleCylinderPortalShell{gctx, vol2};

    BOOST_CHECK_THROW(
        CylinderStackPortalShell(gctx, {&shell1}, BinningValue::binR),
        std::invalid_argument);

    BOOST_CHECK_THROW(
        CylinderStackPortalShell(gctx, {&shell1, &shell2}, BinningValue::binZ),
        std::invalid_argument);
  }
}

// @TODO: Should CylinderStackPortalShell also do the fusing?

BOOST_AUTO_TEST_SUITE_END()  // CylinderStack
BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
