// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <cstddef>
#include <sstream>

using namespace Acts;

namespace ActsTests {

// The binning geometry of the underlying multi-axis (bin lookups, bin edges,
// neighborhoods, closest points, ...) is exercised in `MultiAxisTests.cpp`.
// These tests focus on the behaviour that `Grid` adds on top of its
// multi-axis: value storage and access (`at`, `atPosition`, `atLocalBins`),
// interpolation, type conversion and the `IGrid` interface. Consistency of the
// different value accessors is cross-checked through `grid.multiAxis()`.

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(grid_test_1d_equidistant) {
  using Point = std::array<double, 1>;
  using indices = std::array<std::size_t, 1>;

  const Axis a(0.0, 4.0, 4u);
  Grid g(Type<double>, a);

  BOOST_CHECK_EQUAL(g.size(), 6u);

  // initialize grid
  for (std::size_t bin = 0; bin < g.size(); ++bin) {
    g.at(bin) = bin;
  }

  // consistency of access
  const Point point{0.7};
  const std::size_t globalBin = g.multiAxis().getGlobalBinFromPoint(point);
  const indices localBins = g.multiAxis().getLocalBinsFromGlobalBin(globalBin);

  BOOST_CHECK_EQUAL(g.atPosition(point), g.at(globalBin));
  BOOST_CHECK_EQUAL(g.atPosition(point), g.atLocalBins(localBins));
}

BOOST_AUTO_TEST_CASE(grid_test_2d_equidistant) {
  using Point = std::array<double, 2>;
  using indices = std::array<std::size_t, 2>;

  const Axis a(0.0, 4.0, 4u);
  const Axis b(0.0, 3.0, 3u);
  Grid g(Type<double>, a, b);

  BOOST_CHECK_EQUAL(g.size(), 30u);

  // initialize grid
  for (std::size_t bin = 0; bin < g.size(); ++bin) {
    g.at(bin) = bin;
  }

  // consistency of access
  const Point point{0.7, 1.3};
  const std::size_t globalBin = g.multiAxis().getGlobalBinFromPoint(point);
  const indices localBins = g.multiAxis().getLocalBinsFromGlobalBin(globalBin);

  BOOST_CHECK_EQUAL(g.atPosition(point), g.at(globalBin));
  BOOST_CHECK_EQUAL(g.atPosition(point), g.atLocalBins(localBins));
}

BOOST_AUTO_TEST_CASE(grid_test_3d_equidistant) {
  using Point = std::array<double, 3>;
  using indices = std::array<std::size_t, 3>;

  const Axis a(0.0, 2.0, 2u);
  const Axis b(0.0, 3.0, 3u);
  const Axis c(0.0, 2.0, 2u);
  Grid g(Type<double>, a, b, c);

  BOOST_CHECK_EQUAL(g.size(), 80u);

  // initialize grid
  for (std::size_t bin = 0; bin < g.size(); ++bin) {
    g.at(bin) = bin;
  }

  // consistency of access
  const Point point{0.7, 2.3, 1.3};
  const std::size_t globalBin = g.multiAxis().getGlobalBinFromPoint(point);
  const indices localBins = g.multiAxis().getLocalBinsFromGlobalBin(globalBin);

  BOOST_CHECK_EQUAL(g.atPosition(point), g.at(globalBin));
  BOOST_CHECK_EQUAL(g.atPosition(point), g.atLocalBins(localBins));
}

BOOST_AUTO_TEST_CASE(grid_test_1d_variable) {
  using Point = std::array<double, 1>;
  using indices = std::array<std::size_t, 1>;

  const Axis a({0.0, 1.0, 4.0});
  Grid g(Type<double>, a);

  BOOST_CHECK_EQUAL(g.size(), 4u);

  // initialize grid
  for (std::size_t bin = 0; bin < g.size(); ++bin) {
    g.at(bin) = bin;
  }

  // consistency of access
  const Point point{0.7};
  const std::size_t globalBin = g.multiAxis().getGlobalBinFromPoint(point);
  const indices localBins = g.multiAxis().getLocalBinsFromGlobalBin(globalBin);

  BOOST_CHECK_EQUAL(g.atPosition(point), g.at(globalBin));
  BOOST_CHECK_EQUAL(g.atPosition(point), g.atLocalBins(localBins));
}

BOOST_AUTO_TEST_CASE(grid_test_2d_variable) {
  using Point = std::array<double, 2>;
  using indices = std::array<std::size_t, 2>;

  const Axis a({0.0, 0.5, 3.0});
  const Axis b({0.0, 1.0, 4.0});
  Grid g(Type<double>, a, b);

  BOOST_CHECK_EQUAL(g.size(), 16u);

  // initialize grid
  for (std::size_t bin = 0; bin < g.size(); ++bin) {
    g.at(bin) = bin;
  }

  // consistency of access
  const Point point{0.7, 1.3};
  const std::size_t globalBin = g.multiAxis().getGlobalBinFromPoint(point);
  const indices localBins = g.multiAxis().getLocalBinsFromGlobalBin(globalBin);

  BOOST_CHECK_EQUAL(g.atPosition(point), g.at(globalBin));
  BOOST_CHECK_EQUAL(g.atPosition(point), g.atLocalBins(localBins));
}

BOOST_AUTO_TEST_CASE(grid_test_3d_variable) {
  using Point = std::array<double, 3>;
  using indices = std::array<std::size_t, 3>;

  const Axis a({0.0, 1.0});
  const Axis b({0.0, 0.5, 3.0});
  const Axis c({0.0, 0.5, 3.0, 3.3});
  Grid g(Type<double>, a, b, c);

  BOOST_CHECK_EQUAL(g.size(), 60u);

  // initialize grid
  for (std::size_t bin = 0; bin < g.size(); ++bin) {
    g.at(bin) = bin;
  }

  // consistency of access
  const Point point{0.7, 1.3, 3.7};
  const std::size_t globalBin = g.multiAxis().getGlobalBinFromPoint(point);
  const indices localBins = g.multiAxis().getLocalBinsFromGlobalBin(globalBin);

  BOOST_CHECK_EQUAL(g.atPosition(point), g.at(globalBin));
  BOOST_CHECK_EQUAL(g.atPosition(point), g.atLocalBins(localBins));
}

BOOST_AUTO_TEST_CASE(grid_test_2d_mixed) {
  using Point = std::array<double, 2>;
  using indices = std::array<std::size_t, 2>;

  const Axis a(0.0, 1.0, 4u);
  const Axis b({0.0, 0.5, 3.0});
  Grid g(Type<double>, a, b);

  BOOST_CHECK_EQUAL(g.size(), 24u);

  // initialize grid
  for (std::size_t bin = 0; bin < g.size(); ++bin) {
    g.at(bin) = bin;
  }

  // consistency of access
  const Point point{1.3, 3.7};
  const std::size_t globalBin = g.multiAxis().getGlobalBinFromPoint(point);
  const indices localBins = g.multiAxis().getLocalBinsFromGlobalBin(globalBin);

  BOOST_CHECK_EQUAL(g.atPosition(point), g.at(globalBin));
  BOOST_CHECK_EQUAL(g.atPosition(point), g.atLocalBins(localBins));
}

BOOST_AUTO_TEST_CASE(grid_test_2d_mixed_at) {
  using Point = std::array<double, 2>;

  const Axis a(0.0, 6.0, 4u);
  const Axis b({0.0, 1.5, 3.0});
  Grid g(Type<double>, a, b);

  // initialize the grid
  g.atPosition(Point({{0, 0}})) = 0.;
  g.atPosition(Point({{1.5, 0}})) = 1.;
  g.atPosition(Point({{3, 0}})) = 2.;
  g.atPosition(Point({{4.5, 0}})) = 3.;
  g.atPosition(Point({{6, 0}})) = 4.;
  g.atPosition(Point({{0, 1.5}})) = 5.;
  g.atPosition(Point({{1.5, 1.5}})) = 6.;
  g.atPosition(Point({{3, 1.5}})) = 7.;
  g.atPosition(Point({{4.5, 1.5}})) = 8.;
  g.atPosition(Point({{6, 1.5}})) = 9.;
  g.atPosition(Point({{0, 3}})) = 10.;
  g.atPosition(Point({{1.5, 3}})) = 11.;
  g.atPosition(Point({{3, 3}})) = 12.;
  g.atPosition(Point({{4.5, 3}})) = 13.;
  g.atPosition(Point({{6, 3}})) = 14.;

  // test general properties
  BOOST_CHECK_EQUAL(g.size(), 24u);

  // test some arbitrary points
  BOOST_CHECK_EQUAL(g.atPosition(Point({{1.2, 0.3}})), 0.);
  BOOST_CHECK_EQUAL(g.atPosition(Point({{2.2, 1.3}})), 1.);
  BOOST_CHECK_EQUAL(g.atPosition(Point({{4.9, 1.8}})), 8.);
  BOOST_CHECK_EQUAL(g.atPosition(Point({{3.7, 2.1}})), 7.);
  BOOST_CHECK_EQUAL(g.atPosition(Point({{0.4, 2.3}})), 5.);
}

BOOST_AUTO_TEST_CASE(grid_interpolation) {
  using Point = std::array<double, 3>;

  const Axis a(1.0, 3.0, 2u);
  const Axis b(1.0, 5.0, 2u);
  const Axis c(1.0, 7.0, 2u);
  Grid g(Type<double>, a, b, c);

  g.atPosition(Point{1., 1., 1.}) = 10.;
  g.atPosition(Point{2., 1., 1.}) = 20.;
  g.atPosition(Point{1., 3., 1.}) = 30.;
  g.atPosition(Point{2., 3., 1.}) = 40.;
  g.atPosition(Point{1., 1., 4.}) = 50.;
  g.atPosition(Point{2., 1., 4.}) = 60.;
  g.atPosition(Point{1., 3., 4.}) = 70.;
  g.atPosition(Point{2., 3., 4.}) = 80.;

  CHECK_CLOSE_REL(g.interpolate(Point{1., 1., 1.}), 10., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{2., 1., 1.}), 20., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1., 3., 1.}), 30., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{2., 3., 1.}), 40., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1., 1., 4.}), 50., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{2., 1., 4.}), 60., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1., 3., 4.}), 70., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{2., 3., 4.}), 80., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1.5, 1., 1.}), 15., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1.5, 3., 1.}), 35., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1., 2., 1.}), 20., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{2., 2., 1.}), 30., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1.5, 1., 4.}), 55., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1.5, 3., 4.}), 75., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1., 2., 4.}), 60., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{2., 2., 4.}), 70., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1., 1., 2.5}), 30., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1., 3., 2.5}), 50., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{2., 1., 2.5}), 40., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{2., 3., 2.5}), 60., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1.5, 2., 2.5}), 360. / 8, 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{1.3, 2.1, 1.6}), 32., 1e-6);
  CHECK_CLOSE_REL(g.interpolate(Point{2., 3., 4.}), 80., 1e-6);
}

BOOST_AUTO_TEST_CASE(grid_type_conversion) {
  // Type conversion test
  using Grid2Int =
      Grid<int, Axis<AxisType::Equidistant>, Axis<AxisType::Variable>>;

  const Axis a(0.0, 1.0, 10u);
  const Axis b({0., 1.2, 2.3, 3.4, 4.5, 5.6});
  const Grid g2(Type<double>, a, b);
  const decltype(g2) g2Copy(g2.multiAxis().getAxesTuple());

  static_assert(std::same_as<decltype(g2), decltype(g2Copy)>);

  auto g2ConvertedInt = g2Copy.convertType<int>();
  static_assert(std::same_as<decltype(g2ConvertedInt), Grid2Int>);
}

BOOST_AUTO_TEST_CASE(grid_full_conversion) {
  // The converter class
  struct DoubleToInt {
    // Declare a value type
    using value_type = int;
    // the conversion operator
    int operator()(double d) { return static_cast<int>(d); }
  };

  // Grid conversion test
  const Axis a(0.0, 1.0, 2u);
  Grid g1(Type<double>, a);

  using Point = std::array<double, 1>;
  g1.atPosition(Point{0.3}) = 1.1;
  g1.atPosition(Point{0.6}) = 2.4;

  DoubleToInt d2i;

  auto g1ConvertedInt = g1.convertGrid(d2i);
  BOOST_CHECK_EQUAL(g1ConvertedInt.atPosition(Point{0.3}), 1);
  BOOST_CHECK_EQUAL(g1ConvertedInt.atPosition(Point{0.6}), 2);
}

BOOST_AUTO_TEST_CASE(Output) {
  const Axis a{AxisOpen, 0.0, 1.0, 10u};
  const Axis b{AxisBound, {1, 2, 3}};

  const Grid g(Type<double>, a, b);

  std::stringstream ss;
  ss << g;
  BOOST_CHECK_EQUAL(ss.str(),
                    "Axis<Equidistant, Open>(0, 1, 10, Undefined), "
                    "Axis<Variable, Bound>({1, 2, 3}, Undefined)");

  const IGrid& ig = g;

  ss.str("");

  ss << ig;

  BOOST_CHECK_EQUAL(ss.str(),
                    "Axis<Equidistant, Open>(0, 1, 10, Undefined), "
                    "Axis<Variable, Bound>({1, 2, 3}, Undefined)");
}

BOOST_AUTO_TEST_CASE(Equality) {
  const Axis a{AxisOpen, 0.0, 1.0, 10u};
  const Axis b{AxisBound, {1, 2, 3}};
  const Axis c{AxisClosed, {1, 2, 5}};

  const Grid ab{Type<double>, a, b};
  const Grid ac{Type<double>, a, c};

  BOOST_CHECK_EQUAL(ab, ab);
  BOOST_CHECK_EQUAL(ac, ac);
  BOOST_CHECK_NE(ab, ac);

  const IGrid& iab = ab;
  const IGrid& iac = ac;

  BOOST_CHECK_EQUAL(iab, iab);
  BOOST_CHECK_EQUAL(iac, iac);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
