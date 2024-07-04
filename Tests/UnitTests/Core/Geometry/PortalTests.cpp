// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Surfaces/SurfaceMergingException.hpp"
#include "Acts/Utilities/AxisFwd.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <iostream>
#include <stdexcept>

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

GeometryContext gctx;

BOOST_FIXTURE_TEST_SUITE(Geometry, Fixture)

BOOST_AUTO_TEST_SUITE(GridMerging)
BOOST_AUTO_TEST_SUITE(Merging1dCylinder)

BOOST_AUTO_TEST_CASE(ZDirection) {
  // Need volumes for identity testing, contents should not matter
  auto vol1 = std::make_shared<Volume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(30_mm, 50_mm, 100_mm));
  auto vol2 = std::make_shared<Volume>(
      Transform3::Identity(),
      std::make_shared<CylinderVolumeBounds>(30_mm, 50_mm, 100_mm));

  auto cyl = Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30_mm,
                                                  100_mm);

  BOOST_CHECK_THROW(
      GridPortalLink::make(*cyl, BinningValue::binZ, Axis{AxisBound, 0, 5, 5}),
      std::invalid_argument);

  std::unique_ptr<GridPortalLink> grid1dCyl = GridPortalLink::make(
      *cyl, BinningValue::binZ, Axis{AxisBound, -100_mm, 100_mm, 10});

  // Another cylinder, shifted in z
  auto cyl2 = Surface::makeShared<CylinderSurface>(
      Transform3{Translation3{Vector3::UnitZ() * 150_mm}}, 30_mm, 50_mm);

  std::unique_ptr<GridPortalLink> grid1dCyl2 = GridPortalLink::make(
      *cyl2, BinningValue::binZ, Axis{AxisBound, -50_mm, 50_mm, 5});

  // Completely invalid
  BOOST_CHECK_THROW(
      grid1dCyl->merge(gctx, *grid1dCyl2, BinningValue::binPhi, *logger),
      AssertionFailureException);
  // Invalid direction, as the cylinders are shifted in z, and can't be merged
  // in r x phi
  BOOST_CHECK_THROW(
      grid1dCyl->merge(gctx, *grid1dCyl2, BinningValue::binRPhi, *logger),
      SurfaceMergingException);

  BOOST_TEST_CONTEXT("Consistent equidistant") {
    auto mergedPtr =
        grid1dCyl->merge(gctx, *grid1dCyl2, BinningValue::binZ, *logger);
    BOOST_REQUIRE_NE(mergedPtr, nullptr);
    const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
    BOOST_REQUIRE_NE(merged, nullptr);
    BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
    const auto& axis = *merged->grid().axes().front();
    BOOST_CHECK_EQUAL(axis.getMin(), -150_mm);
    BOOST_CHECK_EQUAL(axis.getMax(), 150_mm);
    BOOST_CHECK_EQUAL(axis.getNBins(), 15);
    BOOST_CHECK_EQUAL(axis.getType(), AxisType::Equidistant);
    BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Bound);
  }

  BOOST_TEST_CONTEXT("Inconsistent equidistant") {
    std::unique_ptr<GridPortalLink> grid1dCyl2BinWidthChanged =
        GridPortalLink::make(*cyl2, BinningValue::binZ,
                             Axis{AxisBound, -50_mm, 50_mm, 6});

    auto mergedPtr = grid1dCyl->merge(gctx, *grid1dCyl2BinWidthChanged,
                                      BinningValue::binZ, *logger);
    BOOST_REQUIRE_NE(mergedPtr, nullptr);
    const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
    BOOST_REQUIRE_NE(merged, nullptr);
    BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
    const auto& axis = *merged->grid().axes().front();
    BOOST_CHECK_EQUAL(axis.getMin(), -150_mm);
    BOOST_CHECK_EQUAL(axis.getMax(), 150_mm);
    BOOST_CHECK_EQUAL(axis.getNBins(), 16);
    BOOST_CHECK_EQUAL(axis.getType(), AxisType::Variable);
    BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Bound);
  }

  BOOST_TEST_CONTEXT("Right Variable") {
    std::unique_ptr<GridPortalLink> gridLeft = GridPortalLink::make(
        *cyl, BinningValue::binZ, Axis{AxisBound, -100_mm, 100_mm, 10});

    std::unique_ptr<GridPortalLink> gridRight =
        GridPortalLink::make(*cyl2, BinningValue::binZ,
                             Axis{AxisBound, {-50_mm, -10_mm, 10_mm, 50_mm}});

    auto mergedPtr =
        gridLeft->merge(gctx, *gridRight, BinningValue::binZ, *logger);
    BOOST_REQUIRE_NE(mergedPtr, nullptr);
    const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
    BOOST_REQUIRE_NE(merged, nullptr);
    BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
    const auto& axis = *merged->grid().axes().front();
    BOOST_CHECK_EQUAL(axis.getMin(), -150_mm);
    BOOST_CHECK_EQUAL(axis.getMax(), 150_mm);
    BOOST_CHECK_EQUAL(axis.getNBins(), 13);
    BOOST_CHECK_EQUAL(axis.getType(), AxisType::Variable);
    BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Bound);
  }

  BOOST_TEST_CONTEXT("Left Variable") {
    std::unique_ptr<GridPortalLink> gridLeft =
        GridPortalLink::make(*cyl, BinningValue::binZ,
                             Axis{AxisBound, {-100_mm, -80_mm, 10_mm, 100_mm}});

    std::unique_ptr<GridPortalLink> gridRight = GridPortalLink::make(
        *cyl2, BinningValue::binZ, Axis{AxisBound, -50_mm, 50_mm, 8});

    auto mergedPtr =
        gridLeft->merge(gctx, *gridRight, BinningValue::binZ, *logger);
    BOOST_REQUIRE_NE(mergedPtr, nullptr);
    const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
    BOOST_REQUIRE_NE(merged, nullptr);
    BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
    const auto& axis = *merged->grid().axes().front();
    BOOST_CHECK_EQUAL(axis.getMin(), -150_mm);
    BOOST_CHECK_EQUAL(axis.getMax(), 150_mm);
    BOOST_CHECK_EQUAL(axis.getNBins(), 11);
    BOOST_CHECK_EQUAL(axis.getType(), AxisType::Variable);
    BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Bound);
  }

  BOOST_TEST_CONTEXT("Both Variable") {
    std::unique_ptr<GridPortalLink> gridLeft =
        GridPortalLink::make(*cyl, BinningValue::binZ,
                             Axis{AxisBound, {-100_mm, -80_mm, 10_mm, 100_mm}});

    std::unique_ptr<GridPortalLink> gridRight =
        GridPortalLink::make(*cyl2, BinningValue::binZ,
                             Axis{AxisBound, {-50_mm, -10_mm, 10_mm, 50_mm}});

    auto mergedPtr =
        gridLeft->merge(gctx, *gridRight, BinningValue::binZ, *logger);
    BOOST_REQUIRE_NE(mergedPtr, nullptr);
    const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
    BOOST_REQUIRE_NE(merged, nullptr);
    BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
    const auto& axis = *merged->grid().axes().front();
    BOOST_CHECK_EQUAL(axis.getMin(), -150_mm);
    BOOST_CHECK_EQUAL(axis.getMax(), 150_mm);
    BOOST_CHECK_EQUAL(axis.getNBins(), 6);
    BOOST_CHECK_EQUAL(axis.getType(), AxisType::Variable);
    BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Bound);
  }

  BOOST_TEST_CONTEXT("Non bound axis") {
    std::unique_ptr<GridPortalLink> gridLeft =
        GridPortalLink::make(*cyl, BinningValue::binZ,
                             Axis{AxisBound, {-100_mm, -80_mm, 10_mm, 100_mm}});
    std::unique_ptr<GridPortalLink> gridRightClosed =
        GridPortalLink::make(*cyl2, BinningValue::binZ,
                             Axis{AxisClosed, {-50_mm, -10_mm, 10_mm, 50_mm}});
    std::unique_ptr<GridPortalLink> gridRightOpen =
        GridPortalLink::make(*cyl2, BinningValue::binZ,
                             Axis{AxisOpen, {-50_mm, -10_mm, 10_mm, 50_mm}});

    // @TODO: Implement fallback to binary for this
    BOOST_CHECK_THROW(
        gridLeft->merge(gctx, *gridRightClosed, BinningValue::binZ, *logger),
        std::logic_error);
    BOOST_CHECK_THROW(
        gridLeft->merge(gctx, *gridRightOpen, BinningValue::binZ, *logger),
        std::logic_error);
  }

  // @TODO: Merge phi sectors with z binning
}

BOOST_AUTO_TEST_CASE(RPhiDirection) {
  auto cyl = Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30_mm,
                                                  100_mm);
  BOOST_CHECK_THROW(GridPortalLink::make(*cyl, BinningValue::binRPhi,
                                         Axis{AxisBound, 0, 5, 5}),
                    std::invalid_argument);

  auto cylNonZeroAverage = Surface::makeShared<CylinderSurface>(
      Transform3::Identity(), 30_mm, 100_mm, 20_degree, 45_degree);

  BOOST_CHECK_THROW(
      GridPortalLink::make(
          *cylNonZeroAverage, BinningValue::binRPhi,
          Axis{AxisBound, -20_degree * 30_mm, 20_degree * 30_mm, 5}),
      std::invalid_argument);

  BOOST_TEST_CONTEXT("Colinear merge in rPhi") {
    auto cylPhi1 = Surface::makeShared<CylinderSurface>(
        Transform3::Identity() * AngleAxis3(45_degree, Vector3::UnitZ()), 30_mm,
        100_mm, 20_degree, 0_degree);

    auto cylPhi2 = Surface::makeShared<CylinderSurface>(
        Transform3::Identity() * AngleAxis3(105_degree, Vector3::UnitZ()),
        30_mm, 100_mm, 40_degree, 0_degree);

    auto portalPhi1 = GridPortalLink::make(
        *cylPhi1, BinningValue::binRPhi,
        Axis{AxisBound, -20_degree * 30_mm, 20_degree * 30_mm, 5});

    auto portalPhi2 = GridPortalLink::make(
        *cylPhi2, BinningValue::binRPhi,
        Axis{AxisBound, -40_degree * 30_mm, 40_degree * 30_mm, 10});

    auto portalMerged =
        portalPhi1->merge(gctx, *portalPhi2, BinningValue::binRPhi, *logger);
    BOOST_REQUIRE_NE(portalMerged, nullptr);

    const auto* merged =
        dynamic_cast<const GridPortalLink*>(portalMerged.get());
    BOOST_REQUIRE_NE(merged, nullptr);
    BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
    const auto& axis = *merged->grid().axes().front();
    BOOST_CHECK_EQUAL(axis.getMin(), -60_degree * 30_mm);
    BOOST_CHECK_EQUAL(axis.getMax(), 60_degree * 30_mm);
    BOOST_CHECK_EQUAL(axis.getNBins(), 15);
    BOOST_CHECK_EQUAL(axis.getType(), AxisType::Equidistant);
    BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Bound);

    // @TODO: Test that if you merge half-circles, we get a closed axis
  }

  // @TODO: Colinear merging of phi sectors along rphi
  // @TODO: Test what happens with surfaces that cover the 0 or 180deg line
  // @TODO: Perpendicular merging (should only be possible with full phi, should not care about axis type)
}
}

BOOST_AUTO_TEST_CASE(Make2DCylinderPortal) {
  auto cyl = Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30_mm,
                                                  100_mm, 45_degree);

  // z bad, rphi bad
  BOOST_CHECK_THROW(GridPortalLink::make(*cyl, Axis{AxisBound, 1, 2, 5},
                                         Axis{AxisBound, 3_mm, 4_mm, 5}),
                    std::invalid_argument);

  // z good, rphi bad
  BOOST_CHECK_THROW(GridPortalLink::make(*cyl, Axis{AxisBound, 3_mm, 4_mm, 5},
                                         Axis{AxisBound, -100_mm, 100_m, 5}),
                    std::invalid_argument);

  // z bad, rphi good
  BOOST_CHECK_THROW(
      GridPortalLink::make(
          *cyl, Axis{AxisBound, -45_degree * 30_mm, 45_degree * 30_mm, 5},
          Axis{AxisBound, -80_mm, 100_mm, 5}),
      std::invalid_argument);

  // z good, rphi good
  BOOST_CHECK_NO_THROW(GridPortalLink::make(
      *cyl, Axis{AxisBound, -45_degree * 30_mm, 45_degree * 30_mm, 5},
      Axis{AxisBound, -100_mm, 100_mm, 5}));
}

// std::unique_ptr<PortalLinkBase> grid1d1 =
//     GridPortalLink::make(Axis{AxisBound{}, 0, 5, 5});
// std::unique_ptr<PortalLinkBase> grid1d2 =
//     GridPortalLink::make(Axis{AxisBound{}, -1, 5, 6});
//
// // @TODO Equidistance different bin widths to variable
// // @TODO Diagonal offset gives error
// // @TODO Zero offset
//
// {
//   auto mergedPtr = grid1d1->merge(*grid1d2, Vector2{5, 0}, *logger);
//   auto* merged = dynamic_cast<GridPortalLink*>(mergedPtr.get());
//   BOOST_REQUIRE_NE(merged, nullptr);
//   BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
//   auto& axis = *merged->grid().axes().front();
//   BOOST_CHECK_EQUAL(axis.getMin(), -1);
//   BOOST_CHECK_EQUAL(axis.getMax(), 10);
// }
//
// {
//   auto mergedPtr = grid1d2->merge(*grid1d1, Vector2{6, 0}, *logger);
//   auto* merged = dynamic_cast<GridPortalLink*>(mergedPtr.get());
//   // BOOST_REQUIRE_NE(merged, nullptr);
//   // merged->grid().axes();
// }

// @TODO Check merge loc0 (closed + bound)
// @TODO Check merge loc1 (bound only)
// @TODO Check non-bound for same direction (error)
// @TODO: Check inconsistent directions (error)

// std::unique_ptr<PortalLinkBase> grid2d1 = GridPortalLink::make(
//     Axis{AxisBound{}, 0, 5, 5}, Axis{AxisBound{}, 0, 5, 5});
//
// grid1d1->merge(*grid2d1);
// grid2d1->merge(*grid1d1);

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
