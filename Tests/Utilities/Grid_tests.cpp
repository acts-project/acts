// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE grid tests
#include <boost/test/included/unit_test.hpp>

#include "ACTS/Utilities/detail/Axis.hpp"
#include "ACTS/Utilities/detail/Grid.hpp"

namespace Acts {

using namespace detail;

namespace Test {

  /// @brief helper class for unit testing Grid
  ///
  /// In order to unit test private methods, a friend helper class is needed.
  template <typename T, class... Axes>
  struct GridTester
  {
    /// @cond
    template <typename... Args>
    GridTester(Args... args) : m_grid(std::forward<Args...>(args...))
    {
    }

    size_t
    getGlobalBinIndex(const std::array<double, sizeof...(Axes)>& point) const
    {
      return m_grid.getGlobalBinIndex(point);
    }

    size_t
    size() const
    {
      return m_grid.size();
    }

  private:
    Grid<T, Axes...> m_grid;
    /// @endcond
  };

  BOOST_AUTO_TEST_CASE(grid_test_1d_equidistant)
  {
    EquidistantAxis a(0.0, 4.0, 4u);
    GridTester<double, EquidistantAxis> g(std::make_tuple(std::move(a)));

    // test general properties
    BOOST_TEST(g.size() == 6u);

    BOOST_TEST(g.getGlobalBinIndex({-0.3}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({-0.}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({0.}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({0.7}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({1}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({1.2}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({2.}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({2.7}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({3.}) == 4u);
    BOOST_TEST(g.getGlobalBinIndex({3.9999}) == 4u);
    BOOST_TEST(g.getGlobalBinIndex({4.}) == 5u);
    BOOST_TEST(g.getGlobalBinIndex({4.98}) == 5u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_2d_equidistant)
  {
    EquidistantAxis a(0.0, 4.0, 4u);
    EquidistantAxis b(0.0, 3.0, 3u);
    GridTester<double, EquidistantAxis, EquidistantAxis> g(
        std::make_tuple(std::move(a), std::move(b)));

    // test general properties
    BOOST_TEST(g.size() == 30u);

    // test grid points
    BOOST_TEST(g.getGlobalBinIndex({-1, -1}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({0, -1}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({1, -1}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({2, -1}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({3, -1}) == 4u);
    BOOST_TEST(g.getGlobalBinIndex({4, -1}) == 5u);
    BOOST_TEST(g.getGlobalBinIndex({-1, 0}) == 6u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0}) == 7u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0}) == 8u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0}) == 9u);
    BOOST_TEST(g.getGlobalBinIndex({3, 0}) == 10u);
    BOOST_TEST(g.getGlobalBinIndex({4, 0}) == 11u);
    BOOST_TEST(g.getGlobalBinIndex({-1, 1}) == 12u);
    BOOST_TEST(g.getGlobalBinIndex({0, 1}) == 13u);
    BOOST_TEST(g.getGlobalBinIndex({1, 1}) == 14u);
    BOOST_TEST(g.getGlobalBinIndex({2, 1}) == 15u);
    BOOST_TEST(g.getGlobalBinIndex({3, 1}) == 16u);
    BOOST_TEST(g.getGlobalBinIndex({4, 1}) == 17u);
    BOOST_TEST(g.getGlobalBinIndex({-1, 2}) == 18u);
    BOOST_TEST(g.getGlobalBinIndex({0, 2}) == 19u);
    BOOST_TEST(g.getGlobalBinIndex({1, 2}) == 20u);
    BOOST_TEST(g.getGlobalBinIndex({2, 2}) == 21u);
    BOOST_TEST(g.getGlobalBinIndex({3, 2}) == 22u);
    BOOST_TEST(g.getGlobalBinIndex({4, 2}) == 23u);
    BOOST_TEST(g.getGlobalBinIndex({-1, 3}) == 24u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3}) == 25u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3}) == 26u);
    BOOST_TEST(g.getGlobalBinIndex({2, 3}) == 27u);
    BOOST_TEST(g.getGlobalBinIndex({3, 3}) == 28u);
    BOOST_TEST(g.getGlobalBinIndex({4, 3}) == 29u);

    // test some arbitrary points
    BOOST_TEST(g.getGlobalBinIndex({1.2, 0.3}) == 8u);
    BOOST_TEST(g.getGlobalBinIndex({2.2, 3.3}) == 27u);
    BOOST_TEST(g.getGlobalBinIndex({0.9, 1.8}) == 13u);
    BOOST_TEST(g.getGlobalBinIndex({3.7, 3.1}) == 28u);
    BOOST_TEST(g.getGlobalBinIndex({1.4, 2.3}) == 20u);
    BOOST_TEST(g.getGlobalBinIndex({-3, 3}) == 24u);
    BOOST_TEST(g.getGlobalBinIndex({8, 1}) == 17u);
    BOOST_TEST(g.getGlobalBinIndex({1, -3}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({3, 11}) == 28u);
    BOOST_TEST(g.getGlobalBinIndex({-2, -3}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({-2, 7}) == 24u);
    BOOST_TEST(g.getGlobalBinIndex({12, -1}) == 5u);
    BOOST_TEST(g.getGlobalBinIndex({12, 11}) == 29u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_3d_equidistant)
  {
    EquidistantAxis a(0.0, 2.0, 2u);
    EquidistantAxis b(0.0, 3.0, 3u);
    EquidistantAxis c(0.0, 2.0, 2u);
    GridTester<double, EquidistantAxis, EquidistantAxis, EquidistantAxis> g(
        std::make_tuple(std::move(a), std::move(b), std::move(c)));

    // test general properties
    BOOST_TEST(g.size() == 80u);

    // test grid points
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 0}) == 25u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 0}) == 26u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0, 0}) == 27u);
    BOOST_TEST(g.getGlobalBinIndex({0, 1, 0}) == 29u);
    BOOST_TEST(g.getGlobalBinIndex({1, 1, 0}) == 30u);
    BOOST_TEST(g.getGlobalBinIndex({2, 1, 0}) == 31u);
    BOOST_TEST(g.getGlobalBinIndex({0, 2, 0}) == 33u);
    BOOST_TEST(g.getGlobalBinIndex({1, 2, 0}) == 34u);
    BOOST_TEST(g.getGlobalBinIndex({2, 2, 0}) == 35u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3, 0}) == 37u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3, 0}) == 38u);
    BOOST_TEST(g.getGlobalBinIndex({2, 3, 0}) == 39u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 1}) == 45u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 1}) == 46u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0, 1}) == 47u);
    BOOST_TEST(g.getGlobalBinIndex({0, 1, 1}) == 49u);
    BOOST_TEST(g.getGlobalBinIndex({1, 1, 1}) == 50u);
    BOOST_TEST(g.getGlobalBinIndex({2, 1, 1}) == 51u);
    BOOST_TEST(g.getGlobalBinIndex({0, 2, 1}) == 53u);
    BOOST_TEST(g.getGlobalBinIndex({1, 2, 1}) == 54u);
    BOOST_TEST(g.getGlobalBinIndex({2, 2, 1}) == 55u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3, 1}) == 57u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3, 1}) == 58u);
    BOOST_TEST(g.getGlobalBinIndex({2, 3, 1}) == 59u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 2}) == 65u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 2}) == 66u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0, 2}) == 67u);
    BOOST_TEST(g.getGlobalBinIndex({0, 1, 2}) == 69u);
    BOOST_TEST(g.getGlobalBinIndex({1, 1, 2}) == 70u);
    BOOST_TEST(g.getGlobalBinIndex({2, 1, 2}) == 71u);
    BOOST_TEST(g.getGlobalBinIndex({0, 2, 2}) == 73u);
    BOOST_TEST(g.getGlobalBinIndex({1, 2, 2}) == 74u);
    BOOST_TEST(g.getGlobalBinIndex({2, 2, 2}) == 75u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3, 2}) == 77u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3, 2}) == 78u);
    BOOST_TEST(g.getGlobalBinIndex({2, 3, 2}) == 79u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_1d_variable)
  {
    VariableAxis a({0.0, 1.0, 4.0});
    GridTester<double, VariableAxis> g(std::make_tuple(std::move(a)));

    // test general properties
    BOOST_TEST(g.size() == 4u);

    BOOST_TEST(g.getGlobalBinIndex({-0.3}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({0.}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({0.7}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({1}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({1.2}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({2.7}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({4.}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({4.98}) == 3u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_2d_variable)
  {
    VariableAxis a({0.0, 1.0, 4.0});
    VariableAxis b({0.0, 0.5, 3.0});
    GridTester<double, VariableAxis, VariableAxis> g(
        std::make_tuple(std::move(a), std::move(b)));

    // test general properties
    BOOST_TEST(g.size() == 16u);

    // test grid points
    BOOST_TEST(g.getGlobalBinIndex({0, 0}) == 5u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0}) == 6u);
    BOOST_TEST(g.getGlobalBinIndex({4, 0}) == 7u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0.5}) == 9u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0.5}) == 10u);
    BOOST_TEST(g.getGlobalBinIndex({4, 0.5}) == 11u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3}) == 13u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3}) == 14u);
    BOOST_TEST(g.getGlobalBinIndex({4, 3}) == 15u);

    // test some arbitrary points
    BOOST_TEST(g.getGlobalBinIndex({1.2, 0.3}) == 6u);
    BOOST_TEST(g.getGlobalBinIndex({2.2, 3.3}) == 14u);
    BOOST_TEST(g.getGlobalBinIndex({0.9, 1.8}) == 9u);
    BOOST_TEST(g.getGlobalBinIndex({0.7, 3.1}) == 13u);
    BOOST_TEST(g.getGlobalBinIndex({1.4, 2.3}) == 10u);
    BOOST_TEST(g.getGlobalBinIndex({-3, 2}) == 8u);
    BOOST_TEST(g.getGlobalBinIndex({8, 1}) == 11u);
    BOOST_TEST(g.getGlobalBinIndex({1, -3}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({3, 11}) == 14u);
    BOOST_TEST(g.getGlobalBinIndex({-2, -3}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({-2, 7}) == 12u);
    BOOST_TEST(g.getGlobalBinIndex({12, -1}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({12, 11}) == 15u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_3d_variable)
  {
    VariableAxis a({0.0, 1.0});
    VariableAxis b({0.0, 0.5, 3.0});
    VariableAxis c({0.0, 0.5, 3.0, 3.3});
    GridTester<double, VariableAxis, VariableAxis, VariableAxis> g(
        std::make_tuple(std::move(a), std::move(b), std::move(c)));

    // test general properties
    BOOST_TEST(g.size() == 60u);

    // test grid points
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 0}) == 16u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 0}) == 17u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0.5, 0}) == 19u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0.5, 0}) == 20u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3, 0}) == 22u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3, 0}) == 23u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 0.5}) == 28u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 0.5}) == 29u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0.5, 0.5}) == 31u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0.5, 0.5}) == 32u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3, 0.5}) == 34u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3, 0.5}) == 35u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 3}) == 40u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 3}) == 41u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0.5, 3}) == 43u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0.5, 3}) == 44u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3, 3}) == 46u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3, 3}) == 47u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 3.3}) == 52u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 3.3}) == 53u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0.5, 3.3}) == 55u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0.5, 3.3}) == 56u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3, 3.3}) == 58u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3, 3.3}) == 59u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_2d_mixed)
  {
    EquidistantAxis a(0.0, 1.0, 4u);
    VariableAxis    b({0.0, 0.5, 3.0});
    GridTester<double, EquidistantAxis, VariableAxis> g(
        std::make_tuple(std::move(a), std::move(b)));

    // test general properties
    BOOST_TEST(g.size() == 24u);

    // test grid points
    BOOST_TEST(g.getGlobalBinIndex({0, 0}) == 7u);
    BOOST_TEST(g.getGlobalBinIndex({0.25, 0}) == 8u);
    BOOST_TEST(g.getGlobalBinIndex({0.5, 0}) == 9u);
    BOOST_TEST(g.getGlobalBinIndex({0.75, 0}) == 10u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0}) == 11u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0.5}) == 13u);
    BOOST_TEST(g.getGlobalBinIndex({0.25, 0.5}) == 14u);
    BOOST_TEST(g.getGlobalBinIndex({0.5, 0.5}) == 15u);
    BOOST_TEST(g.getGlobalBinIndex({0.75, 0.5}) == 16u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0.5}) == 17u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3}) == 19u);
    BOOST_TEST(g.getGlobalBinIndex({0.25, 3}) == 20u);
    BOOST_TEST(g.getGlobalBinIndex({0.5, 3}) == 21u);
    BOOST_TEST(g.getGlobalBinIndex({0.75, 3}) == 22u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3}) == 23u);

    // test some arbitrary points
    BOOST_TEST(g.getGlobalBinIndex({1.2, 0.3}) == 11u);
    BOOST_TEST(g.getGlobalBinIndex({0.2, 1.3}) == 13u);
    BOOST_TEST(g.getGlobalBinIndex({0.9, 1.8}) == 16u);
    BOOST_TEST(g.getGlobalBinIndex({0.7, 2.1}) == 15u);
    BOOST_TEST(g.getGlobalBinIndex({0.4, 0.3}) == 8u);
    BOOST_TEST(g.getGlobalBinIndex({-3, 2}) == 12u);
    BOOST_TEST(g.getGlobalBinIndex({8, 1}) == 17u);
    BOOST_TEST(g.getGlobalBinIndex({0.1, -3}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({0.8, 11}) == 22u);
    BOOST_TEST(g.getGlobalBinIndex({-2, -3}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({-2, 7}) == 18u);
    BOOST_TEST(g.getGlobalBinIndex({12, -1}) == 5u);
    BOOST_TEST(g.getGlobalBinIndex({12, 11}) == 23u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_2d_mixed_at)
  {
    EquidistantAxis a(0.0, 6.0, 4u);
    VariableAxis    b({0.0, 1.5, 3.0});
    Grid<double, EquidistantAxis, VariableAxis> g(
        std::make_tuple(std::move(a), std::move(b)));

    // initialize the grid
    typedef std::array<double, 2> Point;
    g.at(Point({0, 0}))     = 0.;
    g.at(Point({1.5, 0}))   = 1.;
    g.at(Point({3, 0}))     = 2.;
    g.at(Point({4.5, 0}))   = 3.;
    g.at(Point({6, 0}))     = 4.;
    g.at(Point({0, 1.5}))   = 5.;
    g.at(Point({1.5, 1.5})) = 6.;
    g.at(Point({3, 1.5}))   = 7.;
    g.at(Point({4.5, 1.5})) = 8.;
    g.at(Point({6, 1.5}))   = 9.;
    g.at(Point({0, 3}))     = 10.;
    g.at(Point({1.5, 3}))   = 11.;
    g.at(Point({3, 3}))     = 12.;
    g.at(Point({4.5, 3}))   = 13.;
    g.at(Point({6, 3}))     = 14.;

    // test general properties
    BOOST_TEST(g.size() == 24u);

    // test some arbitrary points
    BOOST_TEST(g.at(Point({1.2, 0.3})) == 0.);
    BOOST_TEST(g.at(Point({2.2, 1.3})) == 1.);
    BOOST_TEST(g.at(Point({4.9, 1.8})) == 8.);
    BOOST_TEST(g.at(Point({3.7, 2.1})) == 7.);
    BOOST_TEST(g.at(Point({0.4, 2.3})) == 5.);
  }

}  // namespace Test

}  // namespace Acts
