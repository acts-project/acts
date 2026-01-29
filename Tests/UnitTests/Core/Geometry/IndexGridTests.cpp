// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/IndexGrid.hpp"
#include "Acts/Geometry/ReferenceGenerators.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"

using namespace Acts;

namespace {

/// Helper method to count how many bins are not empty
template <typename indexed_surface_grid>
std::size_t countBins(const indexed_surface_grid& isGrid) {
  std::size_t nonEmptyBins = 0u;
  for (std::size_t igb = 0u; igb < isGrid.grid.size(); ++igb) {
    const auto& gb = isGrid.grid.at(igb);
    if (!gb.empty()) {
      ++nonEmptyBins;
    }
  }
  return nonEmptyBins;
}

}  // namespace

namespace ActsTests {

auto tContext = GeometryContext::dangerouslyDefaultConstruct();

Logging::Level logLevel = Logging::VERBOSE;

BOOST_AUTO_TEST_SUITE(NavigationSuite)

BOOST_AUTO_TEST_CASE(IndexGridXYOneSurfaceCenter) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 0", logLevel));
  ACTS_INFO("Testing X-Y grid.");
  ACTS_INFO("Testing one surface with center generator, should lead to 1 bin.");

  // x-y Axes & Grid
  Axis axisX(AxisBound, -5., 5., 5);
  Axis axisY(AxisBound, -5., 5., 5);
  Grid gridXY(Type<std::vector<unsigned int>>, std::move(axisX),
              std::move(axisY));

  // Indexed Surface grid
  IndexGrid<decltype(gridXY)> indexedGridXY(
      std::move(gridXY), {AxisDirection::AxisX, AxisDirection::AxisY});

  // Create a single surface in the center
  auto rBounds = std::make_shared<RectangleBounds>(4., 4.);
  auto pSurface = Surface::makeShared<PlaneSurface>(Transform3::Identity(),
                                                    std::move(rBounds));

  // The Filler instance and a center based generator
  IndexGridFiller filler{{}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::VERBOSE);
  CenterReferenceGenerator generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {pSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridXY, surfaces, generator);

  std::size_t nonEmptyBins = countBins<decltype(indexedGridXY)>(indexedGridXY);
  // Check the correct number of filled bins
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK_EQUAL(nonEmptyBins, 1u);
}

BOOST_AUTO_TEST_CASE(IndexGridXYOneSurfaceBinValue) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 1", logLevel));
  ACTS_INFO("Testing X-Y grid.");
  ACTS_INFO(
      "Testing one surface with bin value generator, should lead to 1 bin.");

  // x-y Axes & Grid
  Axis axisX(AxisBound, -5., 5., 5);
  Axis axisY(AxisBound, -5., 5., 5);
  Grid gridXY(Type<std::vector<unsigned int>>, std::move(axisX),
              std::move(axisY));

  // Indexed Surface grid
  IndexGrid<decltype(gridXY)> indexedGridXY(
      std::move(gridXY), {AxisDirection::AxisX, AxisDirection::AxisY});

  // Create a single surface in the center
  auto rBounds = std::make_shared<RectangleBounds>(4., 4.);
  auto pSurface = Surface::makeShared<PlaneSurface>(Transform3::Identity(),
                                                    std::move(rBounds));

  // The Filler instance and a center based generator
  IndexGridFiller filler{{}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::VERBOSE);

  AxisDirectionReferenceGenerator<AxisDirection::AxisX> generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {pSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridXY, surfaces, generator);

  std::size_t nonEmptyBins = countBins<decltype(indexedGridXY)>(indexedGridXY);
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK_EQUAL(nonEmptyBins, 1u);
}

BOOST_AUTO_TEST_CASE(IndexGridXYOneSurfacePolyhedron) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 2", logLevel));
  ACTS_INFO("Testing X-Y grid.");
  ACTS_INFO(
      "Testing one surface with polyhedron generator without expansion, should "
      "lead to 5 unique bins, 25 total bins filled");

  // x-y Axes & Grid
  Axis axisX(AxisBound, -5., 5., 5);
  Axis axisY(AxisBound, -5., 5., 5);
  Grid gridXY(Type<std::vector<unsigned int>>, std::move(axisX),
              std::move(axisY));

  // Indexed Surface grid
  IndexGrid<decltype(gridXY)> indexedGridXY(
      std::move(gridXY), {AxisDirection::AxisX, AxisDirection::AxisY});

  // Create a single surface in the center
  auto rBounds = std::make_shared<RectangleBounds>(4., 4.);
  auto pSurface = Surface::makeShared<PlaneSurface>(Transform3::Identity(),
                                                    std::move(rBounds));

  // The Filler instance and a center based generator
  IndexGridFiller filler{{0u, 0u}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::DEBUG);

  PolyhedronReferenceGenerator generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {pSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridXY, surfaces, generator);

  std::size_t nonEmptyBins = countBins<decltype(indexedGridXY)>(indexedGridXY);
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK_EQUAL(nonEmptyBins, 25u);
}

BOOST_AUTO_TEST_CASE(IndexGridXYOneSurfacePolyhedronBinExpansion) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 3", logLevel));
  ACTS_INFO("Testing X-Y grid.");
  ACTS_INFO(
      "Testing one surface with polyhedron generator and expansion, should "
      "lead to 5 unique bins, 49 total bins filled");

  // x-y Axes & Grid
  Axis axisX(AxisBound, -9., 9., 9);
  Axis axisY(AxisBound, -9., 9., 9);
  Grid gridXY(Type<std::vector<unsigned int>>, std::move(axisX),
              std::move(axisY));

  // Indexed Surface grid
  IndexGrid<decltype(gridXY)> indexedGridXY(
      std::move(gridXY), {AxisDirection::AxisX, AxisDirection::AxisY});

  // Create a single surface in the center
  auto rBounds = std::make_shared<RectangleBounds>(4., 4.);
  auto pSurface = Surface::makeShared<PlaneSurface>(Transform3::Identity(),
                                                    std::move(rBounds));

  // The Filler instance and a center based generator
  IndexGridFiller filler{{1u, 1u}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::DEBUG);

  PolyhedronReferenceGenerator generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {pSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridXY, surfaces, generator);

  std::size_t nonEmptyBins = countBins<decltype(indexedGridXY)>(indexedGridXY);
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK_EQUAL(nonEmptyBins, 49u);
}

BOOST_AUTO_TEST_CASE(IndexGridZPhiYOneSurfacePolyhedronBinExpansion) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 4", logLevel));
  ACTS_INFO("Testing Phi-Z grid.");
  ACTS_INFO(
      "Testing one surface with polyhedron generator without expansion, should "
      "lead to 5 unique bins, 6 total bins filled");

  // z-phi Axes & Grid
  Axis axisZ(AxisBound, -9., 9., 9);
  Axis axisPhi(AxisClosed, -std::numbers::pi, std::numbers::pi, 36);
  Grid gridZPhi(Type<std::vector<unsigned int>>, std::move(axisZ),
                std::move(axisPhi));

  // Indexed Surface grid
  IndexGrid<decltype(gridZPhi)> indexedGridZPhi(
      std::move(gridZPhi), {AxisDirection::AxisZ, AxisDirection::AxisPhi});

  auto cBounds =
      std::make_shared<CylinderBounds>(10, 2., std::numbers::pi / 30, 0.);
  auto cSurface = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                       std::move(cBounds));

  // The Filler instance and a center based generator
  IndexGridFiller filler{{0u, 0u}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::DEBUG);

  PolyhedronReferenceGenerator generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {cSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridZPhi, surfaces, generator);

  std::size_t nonEmptyBins =
      countBins<decltype(indexedGridZPhi)>(indexedGridZPhi);
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK_EQUAL(nonEmptyBins, 6u);
}

BOOST_AUTO_TEST_CASE(IndexGridZPhiYOneSurfaceMPIPolyhedronBinExpansion) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 4", logLevel));
  ACTS_INFO("Testing Phi-Z grid.");
  ACTS_INFO(
      "Testing one surface at std::numbers::pi jump, with polyhedron "
      "generator");

  // z-phi Axes & Grid
  Axis axisZ(AxisBound, -9., 9., 9);
  Axis axisPhi(AxisClosed, -std::numbers::pi, std::numbers::pi, 36);
  Grid gridZPhi(Type<std::vector<unsigned int>>, std::move(axisZ),
                std::move(axisPhi));

  // Indexed Surface grid
  IndexGrid<decltype(gridZPhi)> indexedGridZPhi(
      std::move(gridZPhi), {AxisDirection::AxisZ, AxisDirection::AxisPhi});

  auto cBounds =
      std::make_shared<CylinderBounds>(10, 2., std::numbers::pi / 10, 0.);
  auto tf =
      AngleAxis3(std::numbers::pi, Vector3::UnitZ()) * Transform3::Identity();
  auto cSurface = Surface::makeShared<CylinderSurface>(tf, std::move(cBounds));

  // The Filler instance and a center based generator
  IndexGridFiller filler{{0u, 0u}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::DEBUG);

  PolyhedronReferenceGenerator generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {cSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridZPhi, surfaces, generator);

  std::size_t nonEmptyBins =
      countBins<decltype(indexedGridZPhi)>(indexedGridZPhi);
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK_EQUAL(nonEmptyBins, 9u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
