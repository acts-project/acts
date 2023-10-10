// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/detail/IndexedGridFiller.hpp"
#include "Acts/Detector/detail/ReferenceGenerators.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationStateUpdators.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <ostream>
#include <set>
#include <utility>
#include <vector>

using namespace Acts;
using namespace Acts::detail;
using namespace Acts::Experimental;
using namespace Acts::Experimental::detail;

GeometryContext tContext;
Logging::Level logLevel = Logging::VERBOSE;

namespace {

/// Helper method to count how many bins are not empty
template <typename indexed_surface_grid>
std::size_t countBins(const indexed_surface_grid& isGrid) {
  std::size_t nonEmptyBins = 0u;
  for (std::size_t igb = 0u; igb < isGrid.grid.size(); ++igb) {
    const auto& gb = isGrid.grid.at(igb);
    if (not gb.empty()) {
      ++nonEmptyBins;
    }
  }
  return nonEmptyBins;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(BinSequence) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Pre-Test", logLevel));
  ACTS_INFO("Testing bin sequence generators.");

  // Test standard bound local bin sequence
  auto seq48e0b10B = binSequence({4u, 8u}, 0u, 10u, AxisBoundaryType::Bound);
  std::vector<std::size_t> reference = {4u, 5u, 6u, 7u, 8u};
  BOOST_CHECK(seq48e0b10B == reference);

  // Test bound local bin sequence with expansion 1u
  auto seq48e1b10B = binSequence({4u, 8u}, 1u, 10u, AxisBoundaryType::Bound);
  reference = {3u, 4u, 5u, 6u, 7u, 8u, 9u};
  BOOST_CHECK(seq48e1b10B == reference);

  // Test bound local bin sequence with expansion 3u - clipped to max bin 10u
  auto seq48e3b10B = binSequence({4u, 8u}, 3u, 10u, AxisBoundaryType::Bound);
  reference = {1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u};
  BOOST_CHECK(seq48e3b10B == reference);

  // Test open bin sequence with overflow filling
  auto seq48e3b10O = binSequence({4u, 8u}, 3u, 10u, AxisBoundaryType::Open);
  reference = {1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u, 9u, 10u, 11u};
  BOOST_CHECK(seq48e3b10O == reference);

  // Test standard closed local bins
  auto seq48e0b10C = binSequence({4u, 8u}, 0u, 20u, AxisBoundaryType::Closed);
  reference = {4u, 5u, 6u, 7u, 8u};
  BOOST_CHECK(seq48e0b10C == reference);

  // Test closed local bins with expansion
  auto seq48e1b10C = binSequence({4u, 8u}, 1u, 20u, AxisBoundaryType::Closed);
  reference = {3u, 4u, 5u, 6u, 7u, 8u, 9u};
  BOOST_CHECK(seq48e1b10C == reference);

  // Test closed local bins with expansion over bin boundary
  auto seq1029e1b20C =
      binSequence({19u, 20u}, 1u, 20u, AxisBoundaryType::Closed);
  reference = {1u, 18u, 19u, 20u};
  BOOST_CHECK(seq1029e1b20C == reference);

  // Test closed local bins with bin boundary jump
  auto seq218e0b20C = binSequence({2u, 18u}, 0u, 20u, AxisBoundaryType::Closed);
  reference = {1u, 2u, 18u, 19u, 20u};
  BOOST_CHECK(seq218e0b20C == reference);

  // Test closed local bins with bin boundary jump with extension
  auto seq218e2b20C = binSequence({2u, 18u}, 2u, 20u, AxisBoundaryType::Closed);
  reference = {1u, 2u, 3u, 4u, 16u, 17u, 18u, 19u, 20u};
  BOOST_CHECK(seq218e2b20C == reference);
}

BOOST_AUTO_TEST_CASE(IndexGridXYOneSurfaceCenter) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 0", logLevel));
  ACTS_INFO("Testing X-Y grid.");
  ACTS_INFO("Testing one surface with center generator, should lead to 1 bin.");

  // x-y Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisX(-5., 5., 5);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisY(-5., 5., 5);
  Grid<std::vector<unsigned int>, decltype(axisX), decltype(axisY)> gridXY(
      {axisX, axisY});

  // Indexed Surface grid
  IndexedSurfacesImpl<decltype(gridXY)> indexedGridXY(std::move(gridXY),
                                                      {binX, binY});

  // Create a single surface in the center
  auto rBounds = std::make_shared<RectangleBounds>(4., 4.);
  auto pSurface = Surface::makeShared<PlaneSurface>(Transform3::Identity(),
                                                    std::move(rBounds));

  // The Filler instance and a center based generator
  IndexedGridFiller filler{{}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::VERBOSE);
  CenterReferenceGenerator generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {pSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridXY, surfaces, generator);

  std::size_t nonEmptyBins = countBins<decltype(indexedGridXY)>(indexedGridXY);
  // Check the correct number of filled bins
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK(nonEmptyBins == 1u);
}

BOOST_AUTO_TEST_CASE(IndexGridXYOneSurfaceBinValue) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 1", logLevel));
  ACTS_INFO("Testing X-Y grid.");
  ACTS_INFO(
      "Testing one surface with bin value generator, should lead to 1 bin.");

  // x-y Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisX(-5., 5., 5);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisY(-5., 5., 5);
  Grid<std::vector<unsigned int>, decltype(axisX), decltype(axisY)> gridXY(
      {axisX, axisY});

  // Indexed Surface grid
  IndexedSurfacesImpl<decltype(gridXY)> indexedGridXY(std::move(gridXY),
                                                      {binX, binY});

  // Create a single surface in the center
  auto rBounds = std::make_shared<RectangleBounds>(4., 4.);
  auto pSurface = Surface::makeShared<PlaneSurface>(Transform3::Identity(),
                                                    std::move(rBounds));

  // The Filler instance and a center based generator
  IndexedGridFiller filler{{}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::VERBOSE);

  BinningValueReferenceGenerator<binX> generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {pSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridXY, surfaces, generator);

  std::size_t nonEmptyBins = countBins<decltype(indexedGridXY)>(indexedGridXY);
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK(nonEmptyBins == 1u);
}

BOOST_AUTO_TEST_CASE(IndexGridXYOneSurfacePolyhedron) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 2", logLevel));
  ACTS_INFO("Testing X-Y grid.");
  ACTS_INFO(
      "Testing one surface with polyhedron generator without expansion, should "
      "lead to 5 unique bins, 25 total bins filled");

  // x-y Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisX(-5., 5., 5);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisY(-5., 5., 5);
  Grid<std::vector<unsigned int>, decltype(axisX), decltype(axisY)> gridXY(
      {axisX, axisY});

  // Indexed Surface grid
  IndexedSurfacesImpl<decltype(gridXY)> indexedGridXY(std::move(gridXY),
                                                      {binX, binY});

  // Create a single surface in the center
  auto rBounds = std::make_shared<RectangleBounds>(4., 4.);
  auto pSurface = Surface::makeShared<PlaneSurface>(Transform3::Identity(),
                                                    std::move(rBounds));

  // The Filler instance and a center based generator
  IndexedGridFiller filler{{0u, 0u}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::DEBUG);

  PolyhedronReferenceGenerator<1u, true> generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {pSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridXY, surfaces, generator);

  std::size_t nonEmptyBins = countBins<decltype(indexedGridXY)>(indexedGridXY);
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK(nonEmptyBins == 25u);
}

BOOST_AUTO_TEST_CASE(IndexGridXYOneSurfacePolyhedronBinExpansion) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 3", logLevel));
  ACTS_INFO("Testing X-Y grid.");
  ACTS_INFO(
      "Testing one surface with polyhedron generator and expansion, should "
      "lead to 5 unique bins, 49 total bins filled");

  // x-y Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisX(-9., 9., 9);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisY(-9., 9., 9);
  Grid<std::vector<unsigned int>, decltype(axisX), decltype(axisY)> gridXY(
      {axisX, axisY});

  // Indexed Surface grid
  IndexedSurfacesImpl<decltype(gridXY)> indexedGridXY(std::move(gridXY),
                                                      {binX, binY});

  // Create a single surface in the center
  auto rBounds = std::make_shared<RectangleBounds>(4., 4.);
  auto pSurface = Surface::makeShared<PlaneSurface>(Transform3::Identity(),
                                                    std::move(rBounds));

  // The Filler instance and a center based generator
  IndexedGridFiller filler{{1u, 1u}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::DEBUG);

  PolyhedronReferenceGenerator<1u, true> generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {pSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridXY, surfaces, generator);

  std::size_t nonEmptyBins = countBins<decltype(indexedGridXY)>(indexedGridXY);
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK(nonEmptyBins == 49u);
}

BOOST_AUTO_TEST_CASE(IndexGridZPhiYOneSurfacePolyhedronBinExpansion) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 4", logLevel));
  ACTS_INFO("Testing Phi-Z grid.");
  ACTS_INFO(
      "Testing one surface with polyhedron generator without expansion, should "
      "lead to 5 unique bins, 6 total bins filled");

  // z-phi Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisZ(-9., 9., 9);
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(-M_PI, M_PI,
                                                                36);
  Grid<std::vector<unsigned int>, decltype(axisZ), decltype(axisPhi)> gridZPhi(
      {axisZ, axisPhi});

  // Indexed Surface grid
  IndexedSurfacesImpl<decltype(gridZPhi)> indexedGridZPhi(std::move(gridZPhi),
                                                          {binZ, binPhi});

  auto cBounds = std::make_shared<CylinderBounds>(10, 2., M_PI / 30, 0.);
  auto cSurface = Surface::makeShared<CylinderSurface>(Transform3::Identity(),
                                                       std::move(cBounds));

  // The Filler instance and a center based generator
  IndexedGridFiller filler{{0u, 0u}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::DEBUG);

  PolyhedronReferenceGenerator<1u, true> generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {cSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridZPhi, surfaces, generator);

  std::size_t nonEmptyBins =
      countBins<decltype(indexedGridZPhi)>(indexedGridZPhi);
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK(nonEmptyBins == 6u);
}

BOOST_AUTO_TEST_CASE(IndexGridZPhiYOneSurfaceMPIPolyhedronBinExpansion) {
  ACTS_LOCAL_LOGGER(getDefaultLogger("*** Test 4", logLevel));
  ACTS_INFO("Testing Phi-Z grid.");
  ACTS_INFO("Testing one surface at M_PI jump, with polyhedron generator");

  // z-phi Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisZ(-9., 9., 9);
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(-M_PI, M_PI,
                                                                36);
  Grid<std::vector<unsigned int>, decltype(axisZ), decltype(axisPhi)> gridZPhi(
      {axisZ, axisPhi});

  // Indexed Surface grid
  IndexedSurfacesImpl<decltype(gridZPhi)> indexedGridZPhi(std::move(gridZPhi),
                                                          {binZ, binPhi});

  auto cBounds = std::make_shared<CylinderBounds>(10, 2., M_PI / 10, 0.);
  auto tf = AngleAxis3(M_PI, Vector3::UnitZ()) * Transform3::Identity();
  auto cSurface = Surface::makeShared<CylinderSurface>(tf, std::move(cBounds));

  // The Filler instance and a center based generator
  IndexedGridFiller filler{{0u, 0u}};
  filler.oLogger = getDefaultLogger("IndexGridFiller", Logging::DEBUG);

  PolyhedronReferenceGenerator<1u, true> generator;
  std::vector<std::shared_ptr<Surface>> surfaces = {cSurface};

  // Fill the surface
  filler.fill(tContext, indexedGridZPhi, surfaces, generator);

  std::size_t nonEmptyBins =
      countBins<decltype(indexedGridZPhi)>(indexedGridZPhi);
  ACTS_INFO("- filled " << nonEmptyBins << " bins of the grid.");
  BOOST_CHECK(nonEmptyBins == 9u);
}

BOOST_AUTO_TEST_SUITE_END()
