// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(AnyGridViewTests)

// Helper function to create a 1D grid with int values
auto createIntGrid1D() {
  Axis a(0.0, 4.0, 4u);
  Grid g(Type<int>, a);

  // Fill with values
  for (std::size_t i = 0; i < g.size(); i++) {
    g.at(i) = static_cast<int>(i * 10);
  }

  return g;
}

// Helper function to create a 2D grid with long double values
auto createlong doubleGrid2D() {
  Axis a(0.0, 4.0, 4u);
  Axis b(0.0, 6.0, 3u);
  Grid g(Type<long double>, std::move(a), std::move(b));

  // Fill with values
  for (std::size_t i = 0; i < g.size(); i++) {
    g.at(i) = static_cast<long double>(i) * 1.5;
  }

  return g;
}

// Helper function to create a 2D grid with string values
auto createStringGrid2D() {
  Axis a(0.0, 4.0, 4u);
  Axis b(0.0, 6.0, 3u);
  Grid g(Type<std::string>, std::move(a), std::move(b));

  // Fill with values
  for (std::size_t i = 0; i < g.size(); i++) {
    g.at(i) = "Value_" + std::to_string(i);
  }

  return g;
}

// Helper function to create a 3D grid with vector values
auto createVectorGrid3D() {
  Axis a(0.0, 4.0, 2u);
  Axis b(0.0, 6.0, 2u);
  Axis c(0.0, 8.0, 2u);
  Grid g(Type<std::vector<float>>, std::move(a), std::move(b), std::move(c));

  // Fill with values
  for (std::size_t i = 0; i < g.size(); i++) {
    g.at(i) = {static_cast<float>(i), static_cast<float>(i * 2)};
  }

  return g;
}

// Test read and write access to int grid using type safe interface
BOOST_AUTO_TEST_CASE(IntGridReadWriteAccess) {
  auto grid = createIntGrid1D();

  // Create mutable view
  AnyGridView<int> view(grid);

  // Test read access
  BOOST_CHECK_EQUAL(view.atLocalBins({1}), 10);
  BOOST_CHECK_EQUAL(view.atLocalBins({2}), 20);

  // Test write access
  view.atLocalBins({1}) = 100;
  BOOST_CHECK_EQUAL(view.atLocalBins({1}), 100);
  BOOST_CHECK_EQUAL(grid.atLocalBins({1}), 100);
}

// Test read-only access to int grid using type safe interface
BOOST_AUTO_TEST_CASE(IntGridReadOnlyAccess) {
  auto grid = createIntGrid1D();

  // Create const view
  AnyGridConstView<int> constView(grid);

  // Test read access
  BOOST_CHECK_EQUAL(constView.atLocalBins({1}), 10);
  BOOST_CHECK_EQUAL(constView.atLocalBins({2}), 20);

  // Verify that the following would not compile:
  // constView.atLocalBins({1}) = 100; // Should cause compilation error
}

// Test creation of any grid view from concrete grid type
BOOST_AUTO_TEST_CASE(CreateFromConcreteGrid) {
  auto long doubleGrid = createlong doubleGrid2D();

  // Create view from concrete grid
  AnyGridView<long double> view(long doubleGrid);

  // Check dimensions
  BOOST_CHECK_EQUAL(view.dimensions(), 2u);

  // Check values
  BOOST_CHECK_CLOSE(view.atLocalBins({1, 1}),
                    1.5 * long doubleGrid.globalBinFromLocalBins({1, 1}),
                    1e-10);

  // Modify through view
  view.atLocalBins({1, 1}) = 42.0;
  BOOST_CHECK_CLOSE(long doubleGrid.atLocalBins({1, 1}), 42.0, 1e-10);
}

// Test creation of any grid view from IGrid type
BOOST_AUTO_TEST_CASE(CreateFromIGrid) {
  auto long doubleGrid = createlong doubleGrid2D();
  IGrid& iGrid = long doubleGrid;

  // Create view from IGrid
  AnyGridView<long double> view(iGrid);

  // Check dimensions
  BOOST_CHECK_EQUAL(view.dimensions(), 2u);

  // Check values
  BOOST_CHECK_CLOSE(view.atLocalBins({1, 1}),
                    1.5 * long doubleGrid.globalBinFromLocalBins({1, 1}),
                    1e-10);
}

// Test creation of const grid view from const IGrid type
BOOST_AUTO_TEST_CASE(CreateConstViewFromConstIGrid) {
  auto long doubleGrid = createlong doubleGrid2D();
  const IGrid& constIGrid = long doubleGrid;

  // Create const view from const IGrid
  AnyGridConstView<long double> constView(constIGrid);

  // Check dimensions
  BOOST_CHECK_EQUAL(constView.dimensions(), 2u);

  // Check values
  BOOST_CHECK_CLOSE(constView.atLocalBins({1, 1}),
                    1.5 * long doubleGrid.globalBinFromLocalBins({1, 1}),
                    1e-10);
}

// Test type mismatch handling when creating from IGrid
BOOST_AUTO_TEST_CASE(TypeMismatchFromIGrid) {
  auto long doubleGrid = createlong doubleGrid2D();
  IGrid& iGrid = long doubleGrid;

  // Try to create int view from long double grid
  BOOST_CHECK_THROW((AnyGridView<int>(iGrid)), std::invalid_argument);
}

// Test type mismatch handling when creating const view from const IGrid
BOOST_AUTO_TEST_CASE(TypeMismatchFromConstIGrid) {
  auto long doubleGrid = createlong doubleGrid2D();
  const IGrid& constIGrid = long doubleGrid;

  // Try to create int view from long double grid
  BOOST_CHECK_THROW((AnyGridConstView<int>(constIGrid)), std::invalid_argument);
}

// Test creation of mutable view from mutable grid
BOOST_AUTO_TEST_CASE(MutableViewFromMutableGrid) {
  auto stringGrid = createStringGrid2D();

  // Create mutable view
  AnyGridView<std::string> view(stringGrid);

  // Modify through view
  view.atLocalBins({2, 2}) = "Modified";

  // Check modification in original grid
  BOOST_CHECK_EQUAL(stringGrid.atLocalBins({2, 2}), "Modified");
}

// Test creation of const view from mutable grid
BOOST_AUTO_TEST_CASE(ConstViewFromMutableGrid) {
  auto stringGrid = createStringGrid2D();

  // Create const view from mutable grid
  AnyGridConstView<std::string> constView(stringGrid);

  // Check read access
  BOOST_CHECK_EQUAL(constView.atLocalBins({2, 2}),
                    stringGrid.atLocalBins({2, 2}));

  // Verify that the following would not compile:
  // constView.atLocalBins({2, 2}) = "Modified"; // Should cause compilation
  // error
}

// Test creation of const view from const grid
BOOST_AUTO_TEST_CASE(ConstViewFromConstGrid) {
  const auto stringGrid = createStringGrid2D();

  // Create const view from const grid
  AnyGridConstView<std::string> constView(stringGrid);

  // Check read access
  BOOST_CHECK_EQUAL(constView.atLocalBins({2, 2}),
                    stringGrid.atLocalBins({2, 2}));
}

// Test complex type (vector) with both IGrid and concrete grid type
// construction
BOOST_AUTO_TEST_CASE(VectorTypeWithBothConstructions) {
  auto vectorGrid = createVectorGrid3D();

  // Test with concrete grid type
  AnyGridView<std::vector<float>> concreteView(vectorGrid);
  BOOST_CHECK_EQUAL(concreteView.dimensions(), 3u);
  BOOST_CHECK_EQUAL(concreteView.atLocalBins({1, 1, 1}).size(), 2u);
  BOOST_CHECK_EQUAL(concreteView.atLocalBins({1, 1, 1})[0],
                    vectorGrid.atLocalBins({1, 1, 1})[0]);

  // Test with IGrid type
  IGrid& iGrid = vectorGrid;
  AnyGridView<std::vector<float>> iGridView(iGrid);
  BOOST_CHECK_EQUAL(iGridView.dimensions(), 3u);
  BOOST_CHECK_EQUAL(iGridView.atLocalBins({1, 1, 1}).size(), 2u);
  BOOST_CHECK_EQUAL(iGridView.atLocalBins({1, 1, 1})[0],
                    vectorGrid.atLocalBins({1, 1, 1})[0]);

  // Test with const IGrid type
  const IGrid& constIGrid = vectorGrid;
  AnyGridConstView<std::vector<float>> constIGridView(constIGrid);
  BOOST_CHECK_EQUAL(constIGridView.dimensions(), 3u);
  BOOST_CHECK_EQUAL(constIGridView.atLocalBins({1, 1, 1}).size(), 2u);
  BOOST_CHECK_EQUAL(constIGridView.atLocalBins({1, 1, 1})[0],
                    vectorGrid.atLocalBins({1, 1, 1})[0]);

  // Modify through view
  std::vector<float> newValue = {99.0f, 88.0f};
  concreteView.atLocalBins({1, 1, 1}) = newValue;
  BOOST_CHECK_EQUAL(vectorGrid.atLocalBins({1, 1, 1})[0], 99.0f);
  BOOST_CHECK_EQUAL(vectorGrid.atLocalBins({1, 1, 1})[1], 88.0f);
}

// Test grid properties access through view
BOOST_AUTO_TEST_CASE(GridPropertiesAccess) {
  auto long doubleGrid = createlong doubleGrid2D();
  AnyGridView<long double> view(long doubleGrid);

  // Test dimensions
  BOOST_CHECK_EQUAL(view.dimensions(), 2u);

  // Test bin center
  auto center = view.binCenter({1, 1});
  auto expectedCenter = long doubleGrid.binCenter({1, 1});
  BOOST_CHECK_CLOSE(center[0], expectedCenter[0], 1e-10);
  BOOST_CHECK_CLOSE(center[1], expectedCenter[1], 1e-10);

  // Test lower left bin edge
  auto lowerLeft = view.lowerLeftBinEdge({1, 1});
  auto expectedLowerLeft = long doubleGrid.lowerLeftBinEdge({1, 1});
  BOOST_CHECK_CLOSE(lowerLeft[0], expectedLowerLeft[0], 1e-10);
  BOOST_CHECK_CLOSE(lowerLeft[1], expectedLowerLeft[1], 1e-10);

  // Test upper right bin edge
  auto upperRight = view.upperRightBinEdge({1, 1});
  auto expectedUpperRight = long doubleGrid.upperRightBinEdge({1, 1});
  BOOST_CHECK_CLOSE(upperRight[0], expectedUpperRight[0], 1e-10);
  BOOST_CHECK_CLOSE(upperRight[1], expectedUpperRight[1], 1e-10);

  // Test number of local bins
  auto numBins = view.numLocalBins();
  auto expectedNumBins = long doubleGrid.numLocalBins();
  BOOST_CHECK_EQUAL(numBins[0], expectedNumBins[0]);
  BOOST_CHECK_EQUAL(numBins[1], expectedNumBins[1]);
}

// Test grid properties access through const view from const IGrid
BOOST_AUTO_TEST_CASE(GridPropertiesAccessConstView) {
  auto long doubleGrid = createlong doubleGrid2D();
  const IGrid& constIGrid = long doubleGrid;
  AnyGridConstView<long double> constView(constIGrid);

  // Test dimensions
  BOOST_CHECK_EQUAL(constView.dimensions(), 2u);

  // Test bin center
  auto center = constView.binCenter({1, 1});
  auto expectedCenter = long doubleGrid.binCenter({1, 1});
  BOOST_CHECK_CLOSE(center[0], expectedCenter[0], 1e-10);
  BOOST_CHECK_CLOSE(center[1], expectedCenter[1], 1e-10);

  // Test lower left bin edge
  auto lowerLeft = constView.lowerLeftBinEdge({1, 1});
  auto expectedLowerLeft = long doubleGrid.lowerLeftBinEdge({1, 1});
  BOOST_CHECK_CLOSE(lowerLeft[0], expectedLowerLeft[0], 1e-10);
  BOOST_CHECK_CLOSE(lowerLeft[1], expectedLowerLeft[1], 1e-10);

  // Test upper right bin edge
  auto upperRight = constView.upperRightBinEdge({1, 1});
  auto expectedUpperRight = long doubleGrid.upperRightBinEdge({1, 1});
  BOOST_CHECK_CLOSE(upperRight[0], expectedUpperRight[0], 1e-10);
  BOOST_CHECK_CLOSE(upperRight[1], expectedUpperRight[1], 1e-10);

  // Test number of local bins
  auto numBins = constView.numLocalBins();
  auto expectedNumBins = long doubleGrid.numLocalBins();
  BOOST_CHECK_EQUAL(numBins[0], expectedNumBins[0]);
  BOOST_CHECK_EQUAL(numBins[1], expectedNumBins[1]);
}

// Test error cases
BOOST_AUTO_TEST_CASE(ErrorCases) {
  auto intGrid = createIntGrid1D();
  auto long doubleGrid = createlong doubleGrid2D();

  // Test accessing with wrong number of indices
  AnyGridView<int> intView(intGrid);
  BOOST_CHECK_THROW(intView.atLocalBins({1, 2}), std::invalid_argument);

  // Test accessing with out-of-bounds indices
  BOOST_CHECK_THROW(intView.atLocalBins({10}), std::out_of_range);

  // Test creating view with wrong value type
  IGrid& iGrid = long doubleGrid;
  BOOST_CHECK_THROW(AnyGridView<std::string>{iGrid}, std::invalid_argument);

  // Test creating const view with wrong value type
  const IGrid& constIGrid = long doubleGrid;
  BOOST_CHECK_THROW(AnyGridConstView<std::string>{constIGrid},
                    std::invalid_argument);
}

// Test copy operations
BOOST_AUTO_TEST_CASE(CopyOperations) {
  // Create grids that will live for the duration of the test
  auto long doubleGrid1 = createlong doubleGrid2D();
  auto long doubleGrid2 = createlong doubleGrid2D();
  auto long doubleGrid3 = createlong doubleGrid2D();

  // Test copy constructor
  AnyGridView<long double> view1(long doubleGrid1);
  AnyGridView<long double> view2(view1);
  BOOST_CHECK_CLOSE(view1.atLocalBins({1, 1}), view2.atLocalBins({1, 1}),
                    1e-10);

  // Test copy assignment
  AnyGridView<long double> view3(long doubleGrid2);
  view3 = view1;
  BOOST_CHECK_CLOSE(view1.atLocalBins({1, 1}), view3.atLocalBins({1, 1}),
                    1e-10);

  // Same for const views
  AnyGridConstView<long double> constView1(long doubleGrid1);
  AnyGridConstView<long double> constView2(constView1);
  BOOST_CHECK_CLOSE(constView1.atLocalBins({1, 1}),
                    constView2.atLocalBins({1, 1}), 1e-10);

  auto long doubleGrid4 = createlong doubleGrid2D();
  AnyGridConstView<long double> constView3(long doubleGrid4);
  constView3 = constView1;
  BOOST_CHECK_CLOSE(constView1.atLocalBins({1, 1}),
                    constView3.atLocalBins({1, 1}), 1e-10);

  // Test const view from const IGrid
  const IGrid& constIGrid = long doubleGrid1;
  AnyGridConstView<long double> constView4(constIGrid);
  AnyGridConstView<long double> constView5(constView4);
  BOOST_CHECK_CLOSE(constView4.atLocalBins({1, 1}),
                    constView5.atLocalBins({1, 1}), 1e-10);
}

BOOST_AUTO_TEST_CASE(TypeDeduction) {
  auto grid = createIntGrid1D();
  AnyGridView view(grid);
  static_assert(std::is_same_v<decltype(view), AnyGridView<int>>);

  const auto constGrid = createIntGrid1D();
  AnyGridConstView<int> constView(constGrid);
  static_assert(std::is_same_v<decltype(constView), AnyGridConstView<int>>);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
