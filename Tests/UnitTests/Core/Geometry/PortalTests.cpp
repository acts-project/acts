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

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/SurfaceMergingException.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/AxisFwd.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

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

BOOST_AUTO_TEST_SUITE(GridConstruction)

BOOST_AUTO_TEST_CASE(Cylinder) {
  BOOST_TEST_CONTEXT("1D") {
    auto cyl = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                    30_mm, 100_mm);

    // Incompatible binning
    BOOST_CHECK_THROW(GridPortalLink::make(*cyl, BinningValue::binZ,
                                           Axis{AxisBound, 0, 5, 5}),
                      std::invalid_argument);

    std::unique_ptr<GridPortalLink> grid1dCyl = GridPortalLink::make(
        *cyl, BinningValue::binZ, Axis{AxisBound, -100_mm, 100_mm, 10});

    // Throws because non-closed axis
    BOOST_CHECK_THROW(GridPortalLink::make(*cyl, BinningValue::binRPhi,
                                           Axis{AxisBound, -180_degree * 30_mm,
                                                180_degree * 30_mm, 10}),
                      std::invalid_argument);

    std::unique_ptr<GridPortalLink> grid1dCylRPhi = GridPortalLink::make(
        *cyl, BinningValue::binRPhi,
        Axis{AxisClosed, -180_degree * 30_mm, 180_degree * 30_mm, 10});
    BOOST_REQUIRE_NE(grid1dCylRPhi, nullptr);

    Axis axisExpected{AxisClosed, -180_degree * 30_mm, 180_degree * 30_mm, 10};
    BOOST_CHECK_EQUAL(grid1dCylRPhi->grid().axes().size(), 1);
    const auto& axis = *grid1dCylRPhi->grid().axes().front();
    BOOST_CHECK_EQUAL(axis, axisExpected);

    // Another cylinder, shifted in z
    auto cyl2 = Surface::makeShared<CylinderSurface>(
        Transform3{Translation3{Vector3::UnitZ() * 150_mm}}, 30_mm, 50_mm);

    std::unique_ptr<GridPortalLink> grid1dCyl2 = GridPortalLink::make(
        *cyl2, BinningValue::binZ, Axis{AxisBound, -50_mm, 50_mm, 5});

    // Test exception on cylinder with non-zero average phi
    auto cylNonZeroAverage = Surface::makeShared<CylinderSurface>(
        Transform3::Identity(), 30_mm, 100_mm, 20_degree, 45_degree);
    BOOST_CHECK_THROW(
        GridPortalLink::make(*cylNonZeroAverage, BinningValue::binZ,
                             Axis{AxisBound, -100_mm, 100_mm, 10}),
        std::invalid_argument);

    // Extend to a 2D grid with auto phi binning

    auto grid2dCyl1 = grid1dCyl->make2DGrid();
    BOOST_CHECK_EQUAL(grid2dCyl1->grid().axes().size(), 2);
    BOOST_CHECK_EQUAL(grid2dCyl1->surface().bounds(), cyl->bounds());
    const auto& axis1 = *grid2dCyl1->grid().axes().front();
    const auto& axis2 = *grid2dCyl1->grid().axes().back();

    Axis axis1Expected{AxisClosed, -M_PI * 30_mm, M_PI * 30_mm, 1};
    BOOST_CHECK_EQUAL(axis1, axis1Expected);
    Axis axis2Expected{AxisBound, -100_mm, 100_mm, 10};
    BOOST_CHECK_EQUAL(axis2, axis2Expected);

    // @TODO: Check content of the bins

    auto cylPhi = Surface::makeShared<CylinderSurface>(
        Transform3::Identity(), 30_mm, 100_mm, 45_degree);
    std::unique_ptr<GridPortalLink> grid1dCylPhi = GridPortalLink::make(
        *cylPhi, BinningValue::binZ, Axis{AxisBound, -100_mm, 100_mm, 10});

    // Check that phi sector portal does not accept closed axis
    BOOST_CHECK_THROW(GridPortalLink::make(*cylPhi, BinningValue::binRPhi,
                                           Axis{AxisClosed, -45_degree * 30_mm,
                                                45_degree * 30_mm, 10}),
                      std::invalid_argument);

    auto grid2dCylPhi = grid1dCylPhi->make2DGrid();
    BOOST_CHECK_EQUAL(grid2dCylPhi->grid().axes().size(), 2);
    BOOST_CHECK_EQUAL(grid2dCylPhi->surface().bounds(), cylPhi->bounds());
    const auto& axis1Phi = *grid2dCylPhi->grid().axes().front();
    const auto& axis2Phi = *grid2dCylPhi->grid().axes().back();

    Axis axis1PhiExpected{AxisBound, -45_degree * 30_mm, 45_degree * 30_mm, 1};
    BOOST_CHECK_EQUAL(axis1Phi, axis1PhiExpected);
    Axis axis2PhiExpected{AxisBound, -100_mm, 100_mm, 10};
    BOOST_CHECK_EQUAL(axis2Phi, axis2PhiExpected);

    // @TODO: Check content of the bins
  }

  BOOST_TEST_CONTEXT("2D") {
    auto cyl = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                    30_mm, 100_mm, 45_degree);

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

    auto grid1 = GridPortalLink::make(
        *cyl, Axis{AxisBound, -45_degree * 30_mm, 45_degree * 30_mm, 5},
        Axis{AxisBound, -100_mm, 100_mm, 5});

    auto cylFull = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                        30_mm, 100_mm);

    // Throws because non-closed axis
    BOOST_CHECK_THROW(GridPortalLink::make(*cylFull,
                                           Axis{AxisBound, -180_degree * 30_mm,
                                                180_degree * 30_mm, 5},
                                           Axis{AxisBound, -100_mm, 100_mm, 5}),
                      std::invalid_argument);

    auto gridFull = GridPortalLink::make(
        *cylFull, Axis{AxisClosed, -180_degree * 30_mm, 180_degree * 30_mm, 5},
        Axis{AxisBound, -100_mm, 100_mm, 5});

    BOOST_CHECK_EQUAL(gridFull->grid().axes().size(), 2);
    BOOST_CHECK_EQUAL(gridFull->grid().axes().size(), 2);
    Axis axisFullExpected{AxisClosed, -180_degree * 30_mm, 180_degree * 30_mm,
                          5};
    const auto& axisFull = *gridFull->grid().axes().front();
    BOOST_CHECK_EQUAL(axisFull, axisFullExpected);

    // @TODO: Check is same when asked to make 2D grid
  }
}

BOOST_AUTO_TEST_CASE(Disc) {
  BOOST_TEST_CONTEXT("1D") {
    auto disc1 =
        Surface::makeShared<DiscSurface>(Transform3::Identity(), 30_mm, 100_mm);

    BOOST_CHECK_THROW(GridPortalLink::make(*disc1, BinningValue::binZ,
                                           Axis{AxisBound, 30_mm, 100_mm, 3}),
                      std::invalid_argument);

    // Check exception for full disc and non-closed phi axis
    BOOST_CHECK_THROW(
        GridPortalLink::make(*disc1, BinningValue::binPhi,
                             Axis{AxisBound, -180_degree, 180_degree, 3}),
        std::invalid_argument);

    auto grid1 = GridPortalLink::make(*disc1, BinningValue::binR,
                                      Axis{AxisBound, 30_mm, 100_mm, 3});
    BOOST_REQUIRE_NE(grid1, nullptr);
    BOOST_CHECK_EQUAL(grid1->grid().axes().size(), 1);
    const auto& axis = *grid1->grid().axes().front();
    Axis axis1Expected{AxisBound, 30_mm, 100_mm, 3};
    BOOST_CHECK_EQUAL(axis, axis1Expected);

    auto discPhi = Surface::makeShared<DiscSurface>(Transform3::Identity(),
                                                    30_mm, 100_mm, 45_degree);

    // Check thet disc with phi sector does not accept closed axis
    BOOST_CHECK_THROW(
        GridPortalLink::make(*discPhi, BinningValue::binPhi,
                             Axis{AxisClosed, -45_degree, 45_degree, 3}),
        std::invalid_argument);

    auto gridPhi =
        GridPortalLink::make(*discPhi, BinningValue::binPhi,
                             Axis{AxisBound, -45_degree, 45_degree, 3});
    BOOST_REQUIRE_NE(gridPhi, nullptr);

    // Test exception on disc with non-zero average phi
    auto discNonZeroAverage = Surface::makeShared<DiscSurface>(
        Transform3::Identity(),
        std::make_shared<RadialBounds>(30_mm, 100_mm, 45_degree, 75_degree));
    BOOST_CHECK_THROW(
        GridPortalLink::make(*discNonZeroAverage, BinningValue::binR,
                             Axis{AxisBound, 30_mm, 100_mm, 3}),
        std::invalid_argument);

    BOOST_CHECK_EQUAL(gridPhi->grid().axes().size(), 1);
    const auto& axisPhi = *gridPhi->grid().axes().front();
    Axis axisPhi1Expected{AxisBound, -45_degree, 45_degree, 3};
    BOOST_CHECK_EQUAL(axisPhi, axisPhi1Expected);

    // Test making 2D grids from the 1D ones
    // @TODO: Check content of the bins
    auto grid2d = grid1->make2DGrid();
    BOOST_REQUIRE(grid2d);
    BOOST_CHECK_EQUAL(grid2d->grid().axes().size(), 2);
    const auto& axis1 = *grid2d->grid().axes().front();
    BOOST_CHECK_EQUAL(axis1, axis1Expected);
    const auto& axis2 = *grid2d->grid().axes().back();
    BOOST_CHECK_CLOSE(axis2.getMin(), -180_degree, 1e-6);
    BOOST_CHECK_CLOSE(axis2.getMax(), 180_degree, 1e-6);
    BOOST_CHECK_EQUAL(axis2.getNBins(), 1);
    BOOST_CHECK_EQUAL(axis2.getType(), AxisType::Equidistant);
    BOOST_CHECK_EQUAL(axis2.getBoundaryType(), AxisBoundaryType::Closed);

    auto gridPhiBinnedInR = GridPortalLink::make(
        *discPhi, BinningValue::binR, Axis{AxisBound, 30_mm, 100_mm, 3});
    auto grid2dPhiNonClosed = gridPhiBinnedInR->make2DGrid();
    BOOST_REQUIRE(grid2dPhiNonClosed);
    BOOST_CHECK_EQUAL(grid2dPhiNonClosed->grid().axes().size(), 2);
    Axis gridPhiBinnedInRExpected{AxisBound, 30_mm, 100_mm, 3};
    BOOST_CHECK_EQUAL(*grid2dPhiNonClosed->grid().axes().front(),
                      gridPhiBinnedInRExpected);
    const auto& axisPhiNonClosed = *grid2dPhiNonClosed->grid().axes().back();
    BOOST_CHECK_CLOSE(axisPhiNonClosed.getMin(), -45_degree, 1e-6);
    BOOST_CHECK_CLOSE(axisPhiNonClosed.getMax(), 45_degree, 1e-6);
    BOOST_CHECK_EQUAL(axisPhiNonClosed.getNBins(), 1);
    BOOST_CHECK_EQUAL(axisPhiNonClosed.getType(), AxisType::Equidistant);
    BOOST_CHECK_EQUAL(axisPhiNonClosed.getBoundaryType(),
                      AxisBoundaryType::Bound);

    auto grid2dPhi = gridPhi->make2DGrid();
    BOOST_REQUIRE(grid2dPhi);
    BOOST_CHECK_EQUAL(grid2dPhi->grid().axes().size(), 2);
    Axis axis2dPhiExpected{AxisBound, 30_mm, 100_mm, 1};
    BOOST_CHECK_EQUAL(*grid2dPhi->grid().axes().front(), axis2dPhiExpected);
    BOOST_CHECK_EQUAL(*grid2dPhi->grid().axes().back(), axisPhi1Expected);
  }

  BOOST_TEST_CONTEXT("2D") {
    auto discPhi = Surface::makeShared<DiscSurface>(Transform3::Identity(),
                                                    30_mm, 100_mm, 45_degree);

    Axis rBad{AxisBound, 1, 2, 5};
    Axis rGood{AxisBound, 30_mm, 100_mm, 5};
    Axis phiBad{AxisBound, 1, 2, 5};
    Axis phiGood{AxisBound, -45_degree, 45_degree, 5};

    // r bad, phi bad
    BOOST_CHECK_THROW(GridPortalLink::make(*discPhi, rBad, phiBad),
                      std::invalid_argument);
    // r bad, phi good
    BOOST_CHECK_THROW(GridPortalLink::make(*discPhi, rBad, phiGood),
                      std::invalid_argument);
    // r good, phi bad
    BOOST_CHECK_THROW(GridPortalLink::make(*discPhi, rGood, phiBad),
                      std::invalid_argument);
    // r good phi good
    auto grid = GridPortalLink::make(*discPhi, rGood, phiGood);
    BOOST_REQUIRE_NE(grid, nullptr);
  }
}

BOOST_AUTO_TEST_CASE(Plane) {
  // @TODO: Add plane tests
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GridMerging)

BOOST_AUTO_TEST_SUITE(Merging1dCylinder)

BOOST_AUTO_TEST_SUITE(ZDirection)

BOOST_AUTO_TEST_CASE(ColinearMerge) {
  auto cyl = Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30_mm,
                                                  100_mm);

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
}

BOOST_AUTO_TEST_CASE(PerpendicularMerge) {
  // Merge in z direction with phi sectors
  auto cyl1 = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                   30_mm, 100_mm, 35_degree);

  std::unique_ptr<GridPortalLink> grid1 = GridPortalLink::make(
      *cyl1, BinningValue::binRPhi,

      Axis{AxisBound, -35_degree * 30_mm, 35_degree * 30_mm, 3});
  auto cyl2 = Surface::makeShared<CylinderSurface>(
      Transform3{Translation3{Vector3::UnitZ() * 150_mm}}, 30_mm, 50_mm,
      35_degree);

  std::unique_ptr<GridPortalLink> grid2 = GridPortalLink::make(
      *cyl2, BinningValue::binRPhi,
      Axis{AxisBound, -35_degree * 30_mm, 35_degree * 30_mm, 3});

  auto merged12Ptr = grid1->merge(gctx, *grid2, BinningValue::binZ, *logger);
  BOOST_REQUIRE_NE(merged12Ptr, nullptr);
  auto merged12 = dynamic_cast<const GridPortalLink*>(merged12Ptr.get());
  BOOST_REQUIRE_NE(merged12, nullptr);

  BOOST_REQUIRE_EQUAL(merged12->grid().axes().size(), 2);

  const auto& axis1 = *merged12->grid().axes().front();
  const auto& axis2 = *merged12->grid().axes().back();
  Axis axis1Expected{AxisBound, -35_degree * 30_mm, 35_degree * 30_mm, 3};
  BOOST_CHECK_EQUAL(axis1, axis1Expected);
  Axis axis2Expected{AxisBound, {-150_mm, 50_mm, 150_mm}};
  BOOST_CHECK_EQUAL(axis2, axis2Expected);
}

BOOST_AUTO_TEST_SUITE_END()  // ZDirection

BOOST_AUTO_TEST_SUITE(RPhiDirection)

BOOST_AUTO_TEST_CASE(ColinearMerge) {
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

    auto cylPhi3 = Surface::makeShared<CylinderSurface>(
        Transform3::Identity() * AngleAxis3(45_degree, Vector3::UnitZ()), 30_mm,
        100_mm, 90_degree, 0_degree);

    auto cylPhi4 = Surface::makeShared<CylinderSurface>(
        Transform3::Identity() * AngleAxis3(-135_degree, Vector3::UnitZ()),
        30_mm, 100_mm, 90_degree, 0_degree);

    auto portalPhi3 = GridPortalLink::make(
        *cylPhi3, BinningValue::binRPhi,
        Axis{AxisBound, -90_degree * 30_mm, 90_degree * 30_mm, 5});

    auto portalPhi4 = GridPortalLink::make(
        *cylPhi4, BinningValue::binRPhi,
        Axis{AxisBound, -90_degree * 30_mm, 90_degree * 30_mm, 5});

    BOOST_TEST_CONTEXT("Consistent equidistant") {
      auto portalMerged =
          portalPhi1->merge(gctx, *portalPhi2, BinningValue::binRPhi, *logger);
      BOOST_REQUIRE_NE(portalMerged, nullptr);

      const auto* merged =
          dynamic_cast<const GridPortalLink*>(portalMerged.get());
      BOOST_REQUIRE_NE(merged, nullptr);
      BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
      const auto& axis = *merged->grid().axes().front();
      BOOST_CHECK_CLOSE(axis.getMin(), -60_degree * 30_mm, 1e-9);
      BOOST_CHECK_CLOSE(axis.getMax(), 60_degree * 30_mm, 1e-9);
      BOOST_CHECK_EQUAL(axis.getNBins(), 15);
      BOOST_CHECK_EQUAL(axis.getType(), AxisType::Equidistant);
      BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Bound);

      // Test that if you merge half-circles, we get a closed axis
      auto portalMerged34 =
          portalPhi3->merge(gctx, *portalPhi4, BinningValue::binRPhi, *logger);
      BOOST_REQUIRE_NE(portalMerged34, nullptr);

      const auto* merged34 =
          dynamic_cast<const GridPortalLink*>(portalMerged34.get());
      BOOST_REQUIRE_NE(merged34, nullptr);
      BOOST_CHECK_EQUAL(merged34->grid().axes().size(), 1);
      const auto& axis34 = *merged34->grid().axes().front();
      BOOST_CHECK_CLOSE(axis34.getMin(), -180_degree * 30_mm, 1e-9);
      BOOST_CHECK_CLOSE(axis34.getMax(), 180_degree * 30_mm, 1e-9);
      BOOST_CHECK_EQUAL(axis34.getNBins(), 10);
      BOOST_CHECK_EQUAL(axis34.getType(), AxisType::Equidistant);
      BOOST_CHECK_EQUAL(axis34.getBoundaryType(), AxisBoundaryType::Closed);
    }

    BOOST_TEST_CONTEXT("Inconsistent equidistant") {
      auto portalPhi2Mod = GridPortalLink::make(
          *cylPhi2, BinningValue::binRPhi,
          Axis{AxisBound, -40_degree * 30_mm, 40_degree * 30_mm, 3});

      auto portalMergedMod = portalPhi1->merge(gctx, *portalPhi2Mod,
                                               BinningValue::binRPhi, *logger);
      BOOST_REQUIRE_NE(portalMergedMod, nullptr);

      const auto* merged12 =
          dynamic_cast<const GridPortalLink*>(portalMergedMod.get());
      BOOST_REQUIRE_NE(merged12, nullptr);
      BOOST_CHECK_EQUAL(merged12->grid().axes().size(), 1);
      const auto& axis12 = *merged12->grid().axes().front();
      BOOST_CHECK_CLOSE(axis12.getMin(), -60_degree * 30_mm, 1e-9);
      BOOST_CHECK_CLOSE(axis12.getMax(), 60_degree * 30_mm, 1e-9);
      BOOST_CHECK_EQUAL(axis12.getNBins(), 8);
      BOOST_CHECK_EQUAL(axis12.getType(), AxisType::Variable);
      BOOST_CHECK_EQUAL(axis12.getBoundaryType(), AxisBoundaryType::Bound);

      std::vector<ActsScalar> expected12 = {-31.4159, -17.4533, -3.49066,
                                            10.472,   14.6608,  18.8496,
                                            23.0383,  27.2271,  31.4159};
      CHECK_CLOSE_OR_SMALL(axis12.getBinEdges(), expected12, 1e-4, 10e-10);

      auto portalPhi4Mod = GridPortalLink::make(
          *cylPhi4, BinningValue::binRPhi,
          Axis{AxisBound, -90_degree * 30_mm, 90_degree * 30_mm, 3});

      auto portalMerged34 = portalPhi3->merge(gctx, *portalPhi4Mod,
                                              BinningValue::binRPhi, *logger);
      BOOST_REQUIRE_NE(portalMerged34, nullptr);

      const auto* merged34 =
          dynamic_cast<const GridPortalLink*>(portalMerged34.get());
      BOOST_REQUIRE_NE(merged34, nullptr);
      BOOST_CHECK_EQUAL(merged34->grid().axes().size(), 1);
      const auto& axis34 = *merged34->grid().axes().front();
      BOOST_CHECK_CLOSE(axis34.getMin(), -180_degree * 30_mm, 1e-9);
      BOOST_CHECK_CLOSE(axis34.getMax(), 180_degree * 30_mm, 1e-9);
      BOOST_CHECK_EQUAL(axis34.getNBins(), 8);
      BOOST_CHECK_EQUAL(axis34.getType(), AxisType::Variable);
      BOOST_CHECK_EQUAL(axis34.getBoundaryType(), AxisBoundaryType::Closed);

      std::vector<ActsScalar> expected34 = {-94.2478, -75.3982, -56.5487,
                                            -37.6991, -18.8496, 7.10543e-15,
                                            31.4159,  62.8319,  94.2478};
      CHECK_CLOSE_OR_SMALL(axis34.getBinEdges(), expected34, 1e-4, 10e-10);
    }

    BOOST_TEST_CONTEXT("Left variable") {
      BOOST_TEST_CONTEXT("Non-closed") {
        auto gridLeft =
            GridPortalLink::make(*cylPhi1, BinningValue::binRPhi,
                                 Axis{AxisBound,
                                      {-20_degree * 30_mm, -10_degree * 30_mm,
                                       10_degree * 30_mm, 20_degree * 30_mm}});
        auto gridRight = GridPortalLink::make(
            *cylPhi2, BinningValue::binRPhi,
            Axis{AxisBound, -40_degree * 30_mm, 40_degree * 30_mm, 3});

        auto mergedPtr =
            gridLeft->merge(gctx, *gridRight, BinningValue::binRPhi, *logger);
        BOOST_REQUIRE_NE(mergedPtr, nullptr);

        const auto* merged =
            dynamic_cast<const GridPortalLink*>(mergedPtr.get());
        BOOST_REQUIRE_NE(merged, nullptr);
        const auto& axis = *merged->grid().axes().front();
        BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
        BOOST_CHECK_CLOSE(axis.getMin(), -60_degree * 30_mm, 1e-9);
        BOOST_CHECK_CLOSE(axis.getMax(), 60_degree * 30_mm, 1e-9);
        BOOST_CHECK_EQUAL(axis.getNBins(), 6);
        BOOST_CHECK_EQUAL(axis.getType(), AxisType::Variable);
        BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Bound);

        std::vector<ActsScalar> expected = {
            -31.4159, -17.4533, -3.49066, 10.472, 15.708, 26.1799, 31.4159};
        CHECK_CLOSE_OR_SMALL(axis.getBinEdges(), expected, 1e-4, 10e-10);
      }

      BOOST_TEST_CONTEXT("Closed") {
        auto gridLeft = GridPortalLink::make(
            *cylPhi4, BinningValue::binRPhi,
            Axis{AxisBound,
                 {-90_degree * 30_mm, 25_degree * 30_mm, 90_degree * 30_mm}});

        auto gridRight = GridPortalLink::make(
            *cylPhi3, BinningValue::binRPhi,
            Axis{AxisBound, -90_degree * 30_mm, 90_degree * 30_mm, 3});

        auto mergedPtr =
            gridLeft->merge(gctx, *gridRight, BinningValue::binRPhi, *logger);
        BOOST_REQUIRE_NE(mergedPtr, nullptr);

        const auto* merged =
            dynamic_cast<const GridPortalLink*>(mergedPtr.get());
        BOOST_REQUIRE_NE(merged, nullptr);
        BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
        const auto& axis = *merged->grid().axes().front();
        BOOST_CHECK_CLOSE(axis.getMin(), -180_degree * 30_mm, 1e-9);
        BOOST_CHECK_CLOSE(axis.getMax(), 180_degree * 30_mm, 1e-9);
        BOOST_CHECK_EQUAL(axis.getNBins(), 5);
        BOOST_CHECK_EQUAL(axis.getType(), AxisType::Variable);
        BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Closed);

        std::vector<ActsScalar> expected = {-94.2478, -34.0339, 0,
                                            31.4159,  62.8319,  94.2478};
        CHECK_CLOSE_OR_SMALL(axis.getBinEdges(), expected, 1e-4, 10e-10);
      }
    }

    BOOST_TEST_CONTEXT("Right variable") {
      BOOST_TEST_CONTEXT("Non-closed") {
        auto gridLeft = GridPortalLink::make(
            *cylPhi1, BinningValue::binRPhi,
            Axis{AxisBound, -20_degree * 30_mm, 20_degree * 30_mm, 3});
        auto gridRight =
            GridPortalLink::make(*cylPhi2, BinningValue::binRPhi,
                                 Axis{AxisBound,
                                      {-40_degree * 30_mm, -10_degree * 30_mm,
                                       10_degree * 30_mm, 40_degree * 30_mm}});

        auto mergedPtr =
            gridLeft->merge(gctx, *gridRight, BinningValue::binRPhi, *logger);
        BOOST_REQUIRE_NE(mergedPtr, nullptr);

        const auto* merged =
            dynamic_cast<const GridPortalLink*>(mergedPtr.get());
        BOOST_REQUIRE_NE(merged, nullptr);
        const auto& axis = *merged->grid().axes().front();
        BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
        BOOST_CHECK_CLOSE(axis.getMin(), -60_degree * 30_mm, 1e-9);
        BOOST_CHECK_CLOSE(axis.getMax(), 60_degree * 30_mm, 1e-9);
        BOOST_CHECK_EQUAL(axis.getNBins(), 6);
        BOOST_CHECK_EQUAL(axis.getType(), AxisType::Variable);
        BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Bound);

        std::vector<ActsScalar> expected = {-31.4159, -15.708, -5.23599, 10.472,
                                            17.4533,  24.4346, 31.4159};
        CHECK_CLOSE_OR_SMALL(axis.getBinEdges(), expected, 1e-4, 10e-10);
      }

      BOOST_TEST_CONTEXT("Closed") {
        auto gridLeft = GridPortalLink::make(
            *cylPhi4, BinningValue::binRPhi,
            Axis{AxisBound, -90_degree * 30_mm, 90_degree * 30_mm, 3});

        auto gridRight = GridPortalLink::make(
            *cylPhi3, BinningValue::binRPhi,
            Axis{AxisBound,
                 {-90_degree * 30_mm, 25_degree * 30_mm, 90_degree * 30_mm}});

        auto mergedPtr =
            gridLeft->merge(gctx, *gridRight, BinningValue::binRPhi, *logger);
        BOOST_REQUIRE_NE(mergedPtr, nullptr);

        const auto* merged =
            dynamic_cast<const GridPortalLink*>(mergedPtr.get());
        BOOST_REQUIRE_NE(merged, nullptr);
        BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
        const auto& axis = *merged->grid().axes().front();
        BOOST_CHECK_CLOSE(axis.getMin(), -180_degree * 30_mm, 1e-9);
        BOOST_CHECK_CLOSE(axis.getMax(), 180_degree * 30_mm, 1e-9);
        BOOST_CHECK_EQUAL(axis.getNBins(), 5);
        BOOST_CHECK_EQUAL(axis.getType(), AxisType::Variable);
        BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Closed);

        std::vector<ActsScalar> expected = {-94.2478, -62.8319, -31.4159,
                                            0,        60.2139,  94.2478};
        CHECK_CLOSE_OR_SMALL(axis.getBinEdges(), expected, 1e-4, 10e-10);
      }
    }

    BOOST_TEST_CONTEXT("Both variable") {
      BOOST_TEST_CONTEXT("Non-closed") {
        auto gridLeft =
            GridPortalLink::make(*cylPhi1, BinningValue::binRPhi,
                                 Axis{AxisBound,
                                      {-20_degree * 30_mm, -10_degree * 30_mm,
                                       10_degree * 30_mm, 20_degree * 30_mm}});
        auto gridRight = GridPortalLink::make(
            *cylPhi2, BinningValue::binRPhi,
            Axis{AxisBound,
                 {-40_degree * 30_mm, -5_degree * 30_mm, 40_degree * 30_mm}});

        auto mergedPtr =
            gridLeft->merge(gctx, *gridRight, BinningValue::binRPhi, *logger);
        BOOST_REQUIRE_NE(mergedPtr, nullptr);

        const auto* merged =
            dynamic_cast<const GridPortalLink*>(mergedPtr.get());
        BOOST_REQUIRE_NE(merged, nullptr);
        const auto& axis = *merged->grid().axes().front();
        BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
        BOOST_CHECK_CLOSE(axis.getMin(), -60_degree * 30_mm, 1e-9);
        BOOST_CHECK_CLOSE(axis.getMax(), 60_degree * 30_mm, 1e-9);
        BOOST_CHECK_EQUAL(axis.getNBins(), 5);
        BOOST_CHECK_EQUAL(axis.getType(), AxisType::Variable);
        BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Bound);

        std::vector<ActsScalar> expected = {-31.4159, -13.09,  10.472,
                                            15.708,   26.1799, 31.4159};
        CHECK_CLOSE_OR_SMALL(axis.getBinEdges(), expected, 1e-4, 10e-10);
      }

      BOOST_TEST_CONTEXT("Closed") {
        auto gridLeft = GridPortalLink::make(
            *cylPhi4, BinningValue::binRPhi,
            Axis{AxisBound,
                 {-90_degree * 30_mm, 25_degree * 30_mm, 90_degree * 30_mm}});

        auto gridRight =
            GridPortalLink::make(*cylPhi3, BinningValue::binRPhi,
                                 Axis{AxisBound,
                                      {-90_degree * 30_mm, -10_degree * 30_mm,
                                       10_degree * 30_mm, 90_degree * 30_mm}});

        auto mergedPtr =
            gridLeft->merge(gctx, *gridRight, BinningValue::binRPhi, *logger);
        BOOST_REQUIRE_NE(mergedPtr, nullptr);

        const auto* merged =
            dynamic_cast<const GridPortalLink*>(mergedPtr.get());
        BOOST_REQUIRE_NE(merged, nullptr);
        BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
        const auto& axis = *merged->grid().axes().front();
        BOOST_CHECK_CLOSE(axis.getMin(), -180_degree * 30_mm, 1e-9);
        BOOST_CHECK_CLOSE(axis.getMax(), 180_degree * 30_mm, 1e-9);
        BOOST_CHECK_EQUAL(axis.getNBins(), 5);
        BOOST_CHECK_EQUAL(axis.getType(), AxisType::Variable);
        BOOST_CHECK_EQUAL(axis.getBoundaryType(), AxisBoundaryType::Closed);

        std::vector<ActsScalar> expected = {-94.2478, -34.0339, 0,
                                            41.8879,  52.3599,  94.2478};
        CHECK_CLOSE_OR_SMALL(axis.getBinEdges(), expected, 1e-4, 10e-10);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(PerpendicularMerge) {
  // Merge in phi direction with z binning
  auto cylPhi1 = Surface::makeShared<CylinderSurface>(
      Transform3::Identity() * AngleAxis3(45_degree, Vector3::UnitZ()), 30_mm,
      100_mm, 20_degree, 0_degree);

  auto cylPhi2 = Surface::makeShared<CylinderSurface>(
      Transform3::Identity() * AngleAxis3(85_degree, Vector3::UnitZ()), 30_mm,
      100_mm, 20_degree, 0_degree);

  auto portalPhi1 = GridPortalLink::make(*cylPhi1, BinningValue::binZ,
                                         Axis{AxisBound, -100_mm, 100_mm, 5});

  auto portalPhi2 = GridPortalLink::make(*cylPhi2, BinningValue::binZ,
                                         Axis{AxisBound, -100_mm, 100_mm, 5});

  auto merged12Ptr =
      portalPhi1->merge(gctx, *portalPhi2, BinningValue::binRPhi, *logger);
  BOOST_REQUIRE_NE(merged12Ptr, nullptr);
  auto merged12 = dynamic_cast<const GridPortalLink*>(merged12Ptr.get());
  BOOST_REQUIRE_NE(merged12, nullptr);

  const auto& axis1 = *merged12->grid().axes().front();
  const auto& axis2 = *merged12->grid().axes().back();
  // Phi sectors were same size, should give equidistant binning
  Axis axis1Expected{AxisBound, -40_degree * 30_mm, 40_degree * 30_mm, 2};
  BOOST_CHECK_EQUAL(axis1, axis1Expected);
  Axis axis2Expected{AxisBound, -100_mm, 100_mm, 5};
  BOOST_CHECK_EQUAL(axis2, axis2Expected);
}

BOOST_AUTO_TEST_SUITE_END()  // RPhiDirection

BOOST_AUTO_TEST_SUITE_END()  // Merging1dCylinder

BOOST_AUTO_TEST_SUITE(Merging2dCylinder)

BOOST_AUTO_TEST_CASE(ZDirection) {
  BOOST_TEST_CONTEXT("Phi sector") {
    auto cyl1 = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                     30_mm, 100_mm, 45_degree);

    // z good, rphi good
    auto grid1 = GridPortalLink::make(
        *cyl1, Axis{AxisBound, -45_degree * 30_mm, 45_degree * 30_mm, 5},
        Axis{AxisBound, -100_mm, 100_mm, 5});

    auto trf2 = Transform3{Translation3{Vector3::UnitZ() * 150_mm}};
    auto cyl2 =
        Surface::makeShared<CylinderSurface>(trf2, 30_mm, 50_mm, 45_degree);

    // Second grid portal with compatible phi binning
    auto grid2 = GridPortalLink::make(
        *cyl2, Axis{AxisBound, -45_degree * 30_mm, 45_degree * 30_mm, 5},
        Axis{AxisBound, -50_mm, 50_mm, 5});

    // We're merging in z direction, so the phi binnings need to be the same

    auto mergedPtr = grid1->merge(gctx, *grid2, BinningValue::binZ, *logger);
    const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
    BOOST_REQUIRE_NE(mergedPtr, nullptr);

    const auto& axis1 = *merged->grid().axes().front();
    const auto& axis2 = *merged->grid().axes().back();

    BOOST_CHECK_EQUAL(axis1.getMin(), -45_degree * 30_mm);
    BOOST_CHECK_EQUAL(axis1.getMax(), 45_degree * 30_mm);
    BOOST_CHECK_EQUAL(axis1.getNBins(), 5);
    BOOST_CHECK_EQUAL(axis1.getType(), AxisType::Equidistant);
    BOOST_CHECK_EQUAL(axis1.getBoundaryType(), AxisBoundaryType::Bound);

    BOOST_CHECK_EQUAL(axis2.getMin(), -150_mm);
    BOOST_CHECK_EQUAL(axis2.getMax(), 150_mm);
    BOOST_CHECK_EQUAL(axis2.getNBins(), 10);
    BOOST_CHECK_EQUAL(axis2.getType(), AxisType::Variable);
    BOOST_CHECK_EQUAL(axis2.getBoundaryType(), AxisBoundaryType::Bound);

    auto grid3 = GridPortalLink::make(
        *cyl2, Axis{AxisBound, -45_degree * 30_mm, 45_degree * 30_mm, 3},
        Axis{AxisBound, -50_mm, 50_mm, 5});

    // @TODO: Change after binary merging is in
    BOOST_CHECK_THROW(grid1->merge(gctx, *grid3, BinningValue::binZ, *logger),
                      std::logic_error);
  }

  BOOST_TEST_CONTEXT("Check wraparound for full circle in phi") {
    auto cyl1 = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                     30_mm, 100_mm, 180_degree);

    // z good, rphi good
    auto grid1 = GridPortalLink::make(
        *cyl1, Axis{AxisClosed, -180_degree * 30_mm, 180_degree * 30_mm, 5},
        Axis{AxisBound, -100_mm, 100_mm, 5});

    auto trf2 = Transform3{Translation3{Vector3::UnitZ() * 150_mm}};
    auto cyl2 =
        Surface::makeShared<CylinderSurface>(trf2, 30_mm, 50_mm, 180_degree);

    // Second grid portal with compatible phi binning
    auto grid2 = GridPortalLink::make(
        *cyl2, Axis{AxisClosed, -180_degree * 30_mm, 180_degree * 30_mm, 5},
        Axis{AxisBound, -50_mm, 50_mm, 5});

    auto mergedPtr = grid1->merge(gctx, *grid2, BinningValue::binZ, *logger);
    const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
    BOOST_REQUIRE_NE(mergedPtr, nullptr);

    const auto& axis1 = *merged->grid().axes().front();
    const auto& axis2 = *merged->grid().axes().back();

    Axis axis1Expected{AxisClosed, -180_degree * 30_mm, 180_degree * 30_mm, 5};
    BOOST_CHECK_EQUAL(axis1, axis1Expected);
    Axis axis2Expected{AxisBound,
                       {-150, -110, -70, -30, 10, 50, 70, 90, 110, 130, 150}};
    BOOST_CHECK_EQUAL(axis2, axis2Expected);
  }
}

BOOST_AUTO_TEST_CASE(RPhiDirection) {
  auto cyl1 = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                   30_mm, 100_mm, 45_degree);

  // z good, rphi good
  auto grid1 = GridPortalLink::make(
      *cyl1, Axis{AxisBound, -45_degree * 30_mm, 45_degree * 30_mm, 5},
      Axis{AxisBound, -100_mm, 100_mm, 5});
  BOOST_REQUIRE_NE(grid1, nullptr);

  auto trf2 = Transform3{AngleAxis3{90_degree, Vector3::UnitZ()}};
  auto cyl2 =
      Surface::makeShared<CylinderSurface>(trf2, 30_mm, 100_mm, 45_degree);

  // Second grid portal with compatible phi binning
  auto grid2 = GridPortalLink::make(
      *cyl2, Axis{AxisBound, -45_degree * 30_mm, 45_degree * 30_mm, 5},
      Axis{AxisBound, -100_mm, 100_mm, 5});
  BOOST_REQUIRE_NE(grid2, nullptr);

  // We're merging in z direction, so the phi binnings need to be the same

  auto mergedPtr = grid1->merge(gctx, *grid2, BinningValue::binRPhi, *logger);
  const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
  BOOST_REQUIRE_NE(mergedPtr, nullptr);

  const auto& axis1 = *merged->grid().axes().front();
  const auto& axis2 = *merged->grid().axes().back();
  BOOST_CHECK_CLOSE(axis1.getMin(), -90_degree * 30_mm, 1e-8);
  BOOST_CHECK_CLOSE(axis1.getMax(), 90_degree * 30_mm, 1e-8);
  BOOST_CHECK_EQUAL(axis1.getNBins(), 10);
  BOOST_CHECK_EQUAL(axis1.getType(), AxisType::Equidistant);
  BOOST_CHECK_EQUAL(axis1.getBoundaryType(), AxisBoundaryType::Bound);
  Axis axis2Expected{AxisBound, -100_mm, 100_mm, 5};
  BOOST_CHECK_EQUAL(axis2, axis2Expected);
}

BOOST_AUTO_TEST_SUITE_END()  // Merging2dCylinder

BOOST_AUTO_TEST_SUITE(Merging1dDisc)

BOOST_AUTO_TEST_SUITE(RDirection)

BOOST_AUTO_TEST_CASE(ColinearMerge) {
  // Without phi sector
  auto disc1 =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 30_mm, 100_mm);

  auto grid1 = GridPortalLink::make(*disc1, BinningValue::binR,
                                    Axis{AxisBound, 30_mm, 100_mm, 7});

  auto disc2 =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 100_mm, 150_mm);

  auto grid2 = GridPortalLink::make(*disc2, BinningValue::binR,
                                    Axis{AxisBound, 100_mm, 150_mm, 5});

  auto mergedPtr = grid1->merge(gctx, *grid2, BinningValue::binR, *logger);
  BOOST_REQUIRE(mergedPtr);
  const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
  BOOST_REQUIRE_NE(merged, nullptr);

  BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
  Axis axisExpected{AxisBound, 30_mm, 150_mm, 12};
  BOOST_CHECK_EQUAL(*merged->grid().axes().front(), axisExpected);

  // With phi sector
  auto discPhi1 = Surface::makeShared<DiscSurface>(Transform3::Identity(),
                                                   30_mm, 100_mm, 30_degree);

  auto discPhiGrid1 = GridPortalLink::make(*discPhi1, BinningValue::binR,
                                           Axis{AxisBound, 30_mm, 100_mm, 7});

  auto discPhi2 = Surface::makeShared<DiscSurface>(Transform3::Identity(),
                                                   100_mm, 150_mm, 30_degree);

  auto discPhiGrid2 = GridPortalLink::make(*discPhi2, BinningValue::binR,
                                           Axis{AxisBound, 100_mm, 150_mm, 5});

  auto mergedPhiPtr =
      discPhiGrid1->merge(gctx, *discPhiGrid2, BinningValue::binR, *logger);
  BOOST_REQUIRE(mergedPhiPtr);
  const auto* mergedPhi =
      dynamic_cast<const GridPortalLink*>(mergedPhiPtr.get());
  BOOST_REQUIRE_NE(mergedPhi, nullptr);

  BOOST_CHECK_EQUAL(mergedPhi->grid().axes().size(), 1);
  BOOST_CHECK_EQUAL(*mergedPhi->grid().axes().front(), axisExpected);
}

BOOST_AUTO_TEST_CASE(PerpendicularMerge) {
  auto disc1 =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 30_mm, 100_mm);

  auto grid1 =
      GridPortalLink::make(*disc1, BinningValue::binPhi,
                           Axis{AxisClosed, -180_degree, 180_degree, 5});

  auto disc2 =
      Surface::makeShared<DiscSurface>(Transform3::Identity(), 100_mm, 150_mm);

  auto grid2 =
      GridPortalLink::make(*disc2, BinningValue::binPhi,
                           Axis{AxisClosed, -180_degree, 180_degree, 5});

  auto mergedPtr = grid1->merge(gctx, *grid2, BinningValue::binR, *logger);
  BOOST_REQUIRE(mergedPtr);
  const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
  BOOST_REQUIRE_NE(merged, nullptr);

  BOOST_CHECK_EQUAL(merged->grid().axes().size(), 2);
  const auto& axis1 = *merged->grid().axes().front();
  const auto& axis2 = *merged->grid().axes().back();
  Axis axis1Expected{AxisBound, {30_mm, 100_mm, 150_mm}};
  BOOST_CHECK_EQUAL(axis1, axis1Expected);
  Axis axis2Expected{AxisClosed, -180_degree, 180_degree, 5};
  BOOST_CHECK_EQUAL(axis2, axis2Expected);
}

BOOST_AUTO_TEST_SUITE_END()  // RDirection

BOOST_AUTO_TEST_SUITE(PhiDirection)

BOOST_AUTO_TEST_CASE(ColinearMerge) {
  auto disc1 = Surface::makeShared<DiscSurface>(Transform3::Identity(), 30_mm,
                                                100_mm, 30_degree);

  auto grid1 = GridPortalLink::make(*disc1, BinningValue::binPhi,
                                    Axis{AxisBound, -30_degree, 30_degree, 3});

  auto disc2 = Surface::makeShared<DiscSurface>(
      Transform3{AngleAxis3{90_degree, Vector3::UnitZ()}}, 30_mm, 100_mm,
      60_degree);

  auto grid2 = GridPortalLink::make(*disc2, BinningValue::binPhi,
                                    Axis{AxisBound, -60_degree, 60_degree, 6});

  auto mergedPtr = grid1->merge(gctx, *grid2, BinningValue::binPhi, *logger);
  BOOST_REQUIRE(mergedPtr);
  const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
  BOOST_REQUIRE_NE(merged, nullptr);

  BOOST_CHECK_EQUAL(merged->grid().axes().size(), 1);
  const auto& axis = *merged->grid().axes().front();
  BOOST_CHECK_CLOSE(axis.getMin(), -90_degree, 1e-6);
  BOOST_CHECK_CLOSE(axis.getMax(), 90_degree, 1e-6);
  BOOST_CHECK_EQUAL(axis.getNBins(), 9);

  // Check wrapping

  auto disc1Half = Surface::makeShared<DiscSurface>(
      Transform3{AngleAxis3{15_degree, Vector3::UnitZ()}}, 30_mm, 100_mm,
      90_degree);

  auto grid1Half =
      GridPortalLink::make(*disc1Half, BinningValue::binPhi,
                           Axis{AxisBound, -90_degree, 90_degree, 3});

  auto disc2Half = Surface::makeShared<DiscSurface>(
      Transform3{AngleAxis3{-165_degree, Vector3::UnitZ()}}, 30_mm, 100_mm,
      90_degree);

  auto grid2Half =
      GridPortalLink::make(*disc2Half, BinningValue::binPhi,
                           Axis{AxisBound, -90_degree, 90_degree, 3});

  auto mergedHalfPtr =
      grid1Half->merge(gctx, *grid2Half, BinningValue::binPhi, *logger);
  BOOST_REQUIRE(mergedHalfPtr);
  const auto* mergedHalf =
      dynamic_cast<const GridPortalLink*>(mergedHalfPtr.get());
  BOOST_REQUIRE_NE(mergedHalf, nullptr);

  BOOST_CHECK_EQUAL(mergedHalf->grid().axes().size(), 1);
  Axis axisHalfExpected{AxisClosed, -180_degree, 180_degree, 6};
  BOOST_CHECK_EQUAL(axisHalfExpected, *mergedHalf->grid().axes().front());
}

BOOST_AUTO_TEST_CASE(PerpendicularMerge) {
  auto disc1 = Surface::makeShared<DiscSurface>(Transform3::Identity(), 30_mm,
                                                100_mm, 30_degree);

  auto grid1 = GridPortalLink::make(*disc1, BinningValue::binR,
                                    Axis{AxisBound, 30_mm, 100_mm, 3});

  auto disc2 = Surface::makeShared<DiscSurface>(
      Transform3{AngleAxis3{90_degree, Vector3::UnitZ()}}, 30_mm, 100_mm,
      60_degree);

  auto grid2 = GridPortalLink::make(*disc2, BinningValue::binR,
                                    Axis{AxisBound, 30_mm, 100_mm, 3});

  auto mergedPtr = grid1->merge(gctx, *grid2, BinningValue::binPhi, *logger);
  BOOST_REQUIRE(mergedPtr);
  const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
  BOOST_REQUIRE_NE(merged, nullptr);

  BOOST_CHECK_EQUAL(merged->grid().axes().size(), 2);
  const auto& axis1 = *merged->grid().axes().front();
  const auto& axis2 = *merged->grid().axes().back();

  Axis axis1Expected{AxisBound, 30_mm, 100_mm, 3};
  BOOST_CHECK_EQUAL(axis1, axis1Expected);

  BOOST_CHECK_CLOSE(axis2.getMin(), -90_degree, 1e-6);
  BOOST_CHECK_CLOSE(axis2.getMax(), 90_degree, 1e-6);
  BOOST_CHECK_EQUAL(axis2.getNBins(), 2);
  BOOST_CHECK_CLOSE(axis2.getBinEdges().at(1), 30_degree, 1e-6);
  BOOST_CHECK_EQUAL(axis2.getType(), AxisType::Variable);
  BOOST_CHECK_EQUAL(axis2.getBoundaryType(), AxisBoundaryType::Bound);
}

BOOST_AUTO_TEST_SUITE_END()  // PhiDirection

BOOST_AUTO_TEST_SUITE_END()  // Merging2dDisc

BOOST_AUTO_TEST_SUITE(Merging2dDisc)

BOOST_AUTO_TEST_CASE(RDirection) {
  // Basic, because the perpendicular 1D case already tests this to some degree
  auto discPhi1 = Surface::makeShared<DiscSurface>(Transform3::Identity(),
                                                   30_mm, 100_mm, 30_degree);

  auto discPhiGrid1 =
      GridPortalLink::make(*discPhi1, Axis{AxisBound, 30_mm, 100_mm, 7},
                           Axis{AxisBound, -30_degree, 30_degree, 3});

  auto discPhi2 = Surface::makeShared<DiscSurface>(Transform3::Identity(),
                                                   100_mm, 150_mm, 30_degree);

  auto discPhiGrid2 =
      GridPortalLink::make(*discPhi2, Axis{AxisBound, 100_mm, 150_mm, 5},
                           Axis{AxisBound, -30_degree, 30_degree, 3});

  auto mergedPtr =
      discPhiGrid1->merge(gctx, *discPhiGrid2, BinningValue::binR, *logger);
  BOOST_REQUIRE(mergedPtr);
  const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
  BOOST_REQUIRE_NE(merged, nullptr);
  BOOST_CHECK_EQUAL(merged->grid().axes().size(), 2);
  const auto& axis1 = *merged->grid().axes().front();
  const auto& axis2 = *merged->grid().axes().back();

  BOOST_CHECK_EQUAL(axis1.getMin(), 30_mm);
  BOOST_CHECK_EQUAL(axis1.getMax(), 150_mm);
  BOOST_CHECK_EQUAL(axis1.getNBins(), 12);
  BOOST_CHECK_EQUAL(axis1.getType(), AxisType::Equidistant);
  BOOST_CHECK_EQUAL(axis1.getBoundaryType(), AxisBoundaryType::Bound);
  BOOST_CHECK_EQUAL(axis2.getMin(), -30_degree);
  BOOST_CHECK_EQUAL(axis2.getMax(), 30_degree);
  BOOST_CHECK_EQUAL(axis2.getNBins(), 3);
  BOOST_CHECK_EQUAL(axis2.getType(), AxisType::Equidistant);
}

BOOST_AUTO_TEST_CASE(PhiDirection) {
  // Basic, because the perpendicular 1D case already tests this to some degree
  auto disc1 = Surface::makeShared<DiscSurface>(Transform3::Identity(), 30_mm,
                                                100_mm, 30_degree);

  auto grid1 = GridPortalLink::make(*disc1, Axis{AxisBound, 30_mm, 100_mm, 3},
                                    Axis{AxisBound, -30_degree, 30_degree, 3});

  auto disc2 = Surface::makeShared<DiscSurface>(
      Transform3{AngleAxis3{90_degree, Vector3::UnitZ()}}, 30_mm, 100_mm,
      60_degree);

  auto grid2 = GridPortalLink::make(*disc2, Axis{AxisBound, 30_mm, 100_mm, 3},
                                    Axis{AxisBound, -60_degree, 60_degree, 6});

  auto mergedPtr = grid1->merge(gctx, *grid2, BinningValue::binPhi, *logger);
  BOOST_REQUIRE(mergedPtr);
  const auto* merged = dynamic_cast<const GridPortalLink*>(mergedPtr.get());
  BOOST_REQUIRE_NE(merged, nullptr);

  BOOST_CHECK_EQUAL(merged->grid().axes().size(), 2);
  const auto& axis1 = *merged->grid().axes().front();
  const auto& axis2 = *merged->grid().axes().back();

  BOOST_CHECK_EQUAL(axis1.getMin(), 30_mm);
  BOOST_CHECK_EQUAL(axis1.getMax(), 100_mm);
  BOOST_CHECK_EQUAL(axis1.getNBins(), 3);
  BOOST_CHECK_EQUAL(axis1.getType(), AxisType::Equidistant);
  BOOST_CHECK_EQUAL(axis1.getBoundaryType(), AxisBoundaryType::Bound);
  BOOST_CHECK_CLOSE(axis2.getMin(), -90_degree, 1e-6);
  BOOST_CHECK_CLOSE(axis2.getMax(), 90_degree, 1e-6);
  BOOST_CHECK_EQUAL(axis2.getNBins(), 9);
  BOOST_CHECK_EQUAL(axis2.getType(), AxisType::Equidistant);
  BOOST_CHECK_EQUAL(axis2.getBoundaryType(), AxisBoundaryType::Bound);
}

BOOST_AUTO_TEST_SUITE_END()  // Merging2dDisc

BOOST_AUTO_TEST_SUITE_END()  // GridMerging

BOOST_AUTO_TEST_SUITE_END()  // Geometry
}  // namespace Acts::Test
