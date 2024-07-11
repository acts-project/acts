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
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
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

  BOOST_TEST_CONTEXT("Colinear merge") {
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
      std::unique_ptr<GridPortalLink> gridLeft = GridPortalLink::make(
          *cyl, BinningValue::binZ,
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
      std::unique_ptr<GridPortalLink> gridLeft = GridPortalLink::make(
          *cyl, BinningValue::binZ,
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
      std::unique_ptr<GridPortalLink> gridLeft = GridPortalLink::make(
          *cyl, BinningValue::binZ,
          Axis{AxisBound, {-100_mm, -80_mm, 10_mm, 100_mm}});
      std::unique_ptr<GridPortalLink> gridRightClosed = GridPortalLink::make(
          *cyl2, BinningValue::binZ,
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

  BOOST_TEST_CONTEXT("Perpendicular merge") {
    // @TODO: Merge in z direction with phi sectors
  }
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

    auto cylPhi3 = Surface::makeShared<CylinderSurface>(
        Transform3::Identity() * AngleAxis3(45_degree, Vector3::UnitZ()), 30_mm,
        100_mm, 90_degree, 0_degree);

    auto cylPhi4 = Surface::makeShared<CylinderSurface>(
        Transform3::Identity() * AngleAxis3(-135_degree, Vector3::UnitZ()),
        30_mm, 100_mm, 90_degree, 0_degree);

    auto portalPhi3 = GridPortalLink::make(
        *cylPhi3, BinningValue::binRPhi,
        Axis{AxisBound, -90_degree * 30_mm, 90_degree * 30_mm, 5});

    auto portalPhi4Closed = GridPortalLink::make(
        *cylPhi4, BinningValue::binRPhi,
        Axis{AxisClosed, -90_degree * 30_mm, 90_degree * 30_mm, 5});

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

      // @TODO: Cross check binary merging of different boundary types
      // Merge to binary because they have different axis boundary types
      // BOOST_CHECK_THROW(portalPhi3->merge(gctx, *portalPhi4Closed,
      //                                     BinningValue::binRPhi, *logger),
      //                   std::invalid_argument);

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

    // @TODO: Merge in phi direction with z binning
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Merging2dCylinder)

BOOST_AUTO_TEST_CASE(CompatibleMerging) {
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
  auto grid1 = GridPortalLink::make(
      *cyl, Axis{AxisBound, -45_degree * 30_mm, 45_degree * 30_mm, 5},
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
  BOOST_REQUIRE_NE(mergedPtr, nullptr);

  // @TODO: Check wraparound for full circle in phi
  // @TODO: Merge in phi direction with z binning
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
BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
