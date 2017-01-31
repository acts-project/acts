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
    getGlobalBin(const std::array<double, sizeof...(Axes)>& point) const
    {
      return m_grid.getGlobalBin(point);
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
    BOOST_TEST(g.size() == 5u);

    BOOST_TEST(g.getGlobalBin({-0.3}) == 0u);
    BOOST_TEST(g.getGlobalBin({-0.}) == 0u);
    BOOST_TEST(g.getGlobalBin({0.}) == 0u);
    BOOST_TEST(g.getGlobalBin({0.7}) == 0u);
    BOOST_TEST(g.getGlobalBin({1}) == 1u);
    BOOST_TEST(g.getGlobalBin({1.2}) == 1u);
    BOOST_TEST(g.getGlobalBin({2.}) == 2u);
    BOOST_TEST(g.getGlobalBin({2.7}) == 2u);
    BOOST_TEST(g.getGlobalBin({3.}) == 3u);
    BOOST_TEST(g.getGlobalBin({3.9999}) == 3u);
    BOOST_TEST(g.getGlobalBin({4.}) == 4u);
    BOOST_TEST(g.getGlobalBin({4.98}) == 4u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_2d_equidistant)
  {
    EquidistantAxis a(0.0, 4.0, 4u);
    EquidistantAxis b(0.0, 3.0, 3u);
    GridTester<double, EquidistantAxis, EquidistantAxis> g(
        std::make_tuple(std::move(a), std::move(b)));

    // test general properties
    BOOST_TEST(g.size() == 20u);

    // test grid points
    BOOST_TEST(g.getGlobalBin({0, 0}) == 0u);
    BOOST_TEST(g.getGlobalBin({1, 0}) == 1u);
    BOOST_TEST(g.getGlobalBin({2, 0}) == 2u);
    BOOST_TEST(g.getGlobalBin({3, 0}) == 3u);
    BOOST_TEST(g.getGlobalBin({4, 0}) == 4u);
    BOOST_TEST(g.getGlobalBin({0, 1}) == 5u);
    BOOST_TEST(g.getGlobalBin({1, 1}) == 6u);
    BOOST_TEST(g.getGlobalBin({2, 1}) == 7u);
    BOOST_TEST(g.getGlobalBin({3, 1}) == 8u);
    BOOST_TEST(g.getGlobalBin({4, 1}) == 9u);
    BOOST_TEST(g.getGlobalBin({0, 2}) == 10u);
    BOOST_TEST(g.getGlobalBin({1, 2}) == 11u);
    BOOST_TEST(g.getGlobalBin({2, 2}) == 12u);
    BOOST_TEST(g.getGlobalBin({3, 2}) == 13u);
    BOOST_TEST(g.getGlobalBin({4, 2}) == 14u);
    BOOST_TEST(g.getGlobalBin({0, 3}) == 15u);
    BOOST_TEST(g.getGlobalBin({1, 3}) == 16u);
    BOOST_TEST(g.getGlobalBin({2, 3}) == 17u);
    BOOST_TEST(g.getGlobalBin({3, 3}) == 18u);
    BOOST_TEST(g.getGlobalBin({4, 3}) == 19u);

    // test some arbitrary points
    BOOST_TEST(g.getGlobalBin({1.2, 0.3}) == 1u);
    BOOST_TEST(g.getGlobalBin({2.2, 3.3}) == 17u);
    BOOST_TEST(g.getGlobalBin({0.9, 1.8}) == 5u);
    BOOST_TEST(g.getGlobalBin({3.7, 3.1}) == 18u);
    BOOST_TEST(g.getGlobalBin({1.4, 2.3}) == 11u);

    // some outside points
    BOOST_TEST(g.getGlobalBin({-3, 3}) == 15u);
    BOOST_TEST(g.getGlobalBin({8, 1}) == 9u);
    BOOST_TEST(g.getGlobalBin({1, -3}) == 1u);
    BOOST_TEST(g.getGlobalBin({3, 11}) == 18u);
    BOOST_TEST(g.getGlobalBin({-2, -3}) == 0u);
    BOOST_TEST(g.getGlobalBin({-2, 7}) == 15u);
    BOOST_TEST(g.getGlobalBin({12, -1}) == 4u);
    BOOST_TEST(g.getGlobalBin({12, 11}) == 19u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_3d_equidistant)
  {
    EquidistantAxis a(0.0, 2.0, 2u);
    EquidistantAxis b(0.0, 3.0, 3u);
    EquidistantAxis c(0.0, 2.0, 2u);
    GridTester<double, EquidistantAxis, EquidistantAxis, EquidistantAxis> g(
        std::make_tuple(std::move(a), std::move(b), std::move(c)));

    // test general properties
    BOOST_TEST(g.size() == 36u);

    // test grid points
    BOOST_TEST(g.getGlobalBin({0, 0, 0}) == 0u);
    BOOST_TEST(g.getGlobalBin({1, 0, 0}) == 1u);
    BOOST_TEST(g.getGlobalBin({2, 0, 0}) == 2u);
    BOOST_TEST(g.getGlobalBin({0, 1, 0}) == 3u);
    BOOST_TEST(g.getGlobalBin({1, 1, 0}) == 4u);
    BOOST_TEST(g.getGlobalBin({2, 1, 0}) == 5u);
    BOOST_TEST(g.getGlobalBin({0, 2, 0}) == 6u);
    BOOST_TEST(g.getGlobalBin({1, 2, 0}) == 7u);
    BOOST_TEST(g.getGlobalBin({2, 2, 0}) == 8u);
    BOOST_TEST(g.getGlobalBin({0, 3, 0}) == 9u);
    BOOST_TEST(g.getGlobalBin({1, 3, 0}) == 10u);
    BOOST_TEST(g.getGlobalBin({2, 3, 0}) == 11u);
    BOOST_TEST(g.getGlobalBin({0, 0, 1}) == 12u);
    BOOST_TEST(g.getGlobalBin({1, 0, 1}) == 13u);
    BOOST_TEST(g.getGlobalBin({2, 0, 1}) == 14u);
    BOOST_TEST(g.getGlobalBin({0, 1, 1}) == 15u);
    BOOST_TEST(g.getGlobalBin({1, 1, 1}) == 16u);
    BOOST_TEST(g.getGlobalBin({2, 1, 1}) == 17u);
    BOOST_TEST(g.getGlobalBin({0, 2, 1}) == 18u);
    BOOST_TEST(g.getGlobalBin({1, 2, 1}) == 19u);
    BOOST_TEST(g.getGlobalBin({2, 2, 1}) == 20u);
    BOOST_TEST(g.getGlobalBin({0, 3, 1}) == 21u);
    BOOST_TEST(g.getGlobalBin({1, 3, 1}) == 22u);
    BOOST_TEST(g.getGlobalBin({2, 3, 1}) == 23u);
    BOOST_TEST(g.getGlobalBin({0, 0, 2}) == 24u);
    BOOST_TEST(g.getGlobalBin({1, 0, 2}) == 25u);
    BOOST_TEST(g.getGlobalBin({2, 0, 2}) == 26u);
    BOOST_TEST(g.getGlobalBin({0, 1, 2}) == 27u);
    BOOST_TEST(g.getGlobalBin({1, 1, 2}) == 28u);
    BOOST_TEST(g.getGlobalBin({2, 1, 2}) == 29u);
    BOOST_TEST(g.getGlobalBin({0, 2, 2}) == 30u);
    BOOST_TEST(g.getGlobalBin({1, 2, 2}) == 31u);
    BOOST_TEST(g.getGlobalBin({2, 2, 2}) == 32u);
    BOOST_TEST(g.getGlobalBin({0, 3, 2}) == 33u);
    BOOST_TEST(g.getGlobalBin({1, 3, 2}) == 34u);
    BOOST_TEST(g.getGlobalBin({2, 3, 2}) == 35u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_1d_variable)
  {
    VariableAxis a({0.0, 1.0, 4.0});
    GridTester<double, VariableAxis> g(std::make_tuple(std::move(a)));

    // test general properties
    BOOST_TEST(g.size() == 3u);

    BOOST_TEST(g.getGlobalBin({-0.3}) == 0u);
    BOOST_TEST(g.getGlobalBin({0.}) == 0u);
    BOOST_TEST(g.getGlobalBin({0.7}) == 0u);
    BOOST_TEST(g.getGlobalBin({1}) == 1u);
    BOOST_TEST(g.getGlobalBin({1.2}) == 1u);
    BOOST_TEST(g.getGlobalBin({2.7}) == 1u);
    BOOST_TEST(g.getGlobalBin({4.}) == 2u);
    BOOST_TEST(g.getGlobalBin({4.98}) == 2u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_2d_variable)
  {
    VariableAxis a({0.0, 1.0, 4.0});
    VariableAxis b({0.0, 0.5, 3.0});
    GridTester<double, VariableAxis, VariableAxis> g(
        std::make_tuple(std::move(a), std::move(b)));

    // test general properties
    BOOST_TEST(g.size() == 9u);

    // test grid points
    BOOST_TEST(g.getGlobalBin({0, 0}) == 0u);
    BOOST_TEST(g.getGlobalBin({1, 0}) == 1u);
    BOOST_TEST(g.getGlobalBin({4, 0}) == 2u);
    BOOST_TEST(g.getGlobalBin({0, 0.5}) == 3u);
    BOOST_TEST(g.getGlobalBin({1, 0.5}) == 4u);
    BOOST_TEST(g.getGlobalBin({4, 0.5}) == 5u);
    BOOST_TEST(g.getGlobalBin({0, 3}) == 6u);
    BOOST_TEST(g.getGlobalBin({1, 3}) == 7u);
    BOOST_TEST(g.getGlobalBin({4, 3}) == 8u);

    // test some arbitrary points
    BOOST_TEST(g.getGlobalBin({1.2, 0.3}) == 1u);
    BOOST_TEST(g.getGlobalBin({2.2, 3.3}) == 7u);
    BOOST_TEST(g.getGlobalBin({0.9, 1.8}) == 3u);
    BOOST_TEST(g.getGlobalBin({0.7, 3.1}) == 6u);
    BOOST_TEST(g.getGlobalBin({1.4, 2.3}) == 4u);

    // some outside points
    BOOST_TEST(g.getGlobalBin({-3, 2}) == 3u);
    BOOST_TEST(g.getGlobalBin({8, 1}) == 5u);
    BOOST_TEST(g.getGlobalBin({1, -3}) == 1u);
    BOOST_TEST(g.getGlobalBin({3, 11}) == 7u);
    BOOST_TEST(g.getGlobalBin({-2, -3}) == 0u);
    BOOST_TEST(g.getGlobalBin({-2, 7}) == 6u);
    BOOST_TEST(g.getGlobalBin({12, -1}) == 2u);
    BOOST_TEST(g.getGlobalBin({12, 11}) == 8u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_3d_variable)
  {
    VariableAxis a({0.0, 1.0});
    VariableAxis b({0.0, 0.5, 3.0});
    VariableAxis c({0.0, 0.5, 3.0, 3.3});
    GridTester<double, VariableAxis, VariableAxis, VariableAxis> g(
        std::make_tuple(std::move(a), std::move(b), std::move(c)));

    // test general properties
    BOOST_TEST(g.size() == 24u);

    // test grid points
    BOOST_TEST(g.getGlobalBin({0, 0, 0}) == 0u);
    BOOST_TEST(g.getGlobalBin({1, 0, 0}) == 1u);
    BOOST_TEST(g.getGlobalBin({0, 0.5, 0}) == 2u);
    BOOST_TEST(g.getGlobalBin({1, 0.5, 0}) == 3u);
    BOOST_TEST(g.getGlobalBin({0, 3, 0}) == 4u);
    BOOST_TEST(g.getGlobalBin({1, 3, 0}) == 5u);
    BOOST_TEST(g.getGlobalBin({0, 0, 0.5}) == 6u);
    BOOST_TEST(g.getGlobalBin({1, 0, 0.5}) == 7u);
    BOOST_TEST(g.getGlobalBin({0, 0.5, 0.5}) == 8u);
    BOOST_TEST(g.getGlobalBin({1, 0.5, 0.5}) == 9u);
    BOOST_TEST(g.getGlobalBin({0, 3, 0.5}) == 10u);
    BOOST_TEST(g.getGlobalBin({1, 3, 0.5}) == 11u);
    BOOST_TEST(g.getGlobalBin({0, 0, 3}) == 12u);
    BOOST_TEST(g.getGlobalBin({1, 0, 3}) == 13u);
    BOOST_TEST(g.getGlobalBin({0, 0.5, 3}) == 14u);
    BOOST_TEST(g.getGlobalBin({1, 0.5, 3}) == 15u);
    BOOST_TEST(g.getGlobalBin({0, 3, 3}) == 16u);
    BOOST_TEST(g.getGlobalBin({1, 3, 3}) == 17u);
    BOOST_TEST(g.getGlobalBin({0, 0, 3.3}) == 18u);
    BOOST_TEST(g.getGlobalBin({1, 0, 3.3}) == 19u);
    BOOST_TEST(g.getGlobalBin({0, 0.5, 3.3}) == 20u);
    BOOST_TEST(g.getGlobalBin({1, 0.5, 3.3}) == 21u);
    BOOST_TEST(g.getGlobalBin({0, 3, 3.3}) == 22u);
    BOOST_TEST(g.getGlobalBin({1, 3, 3.3}) == 23u);
  }

  BOOST_AUTO_TEST_CASE(grid_test_2d_mixed)
  {
    EquidistantAxis a(0.0, 1.0, 4u);
    VariableAxis    b({0.0, 0.5, 3.0});
    GridTester<double, EquidistantAxis, VariableAxis> g(
        std::make_tuple(std::move(a), std::move(b)));

    // test general properties
    BOOST_TEST(g.size() == 15u);

    // test grid points
    BOOST_TEST(g.getGlobalBin({0, 0}) == 0u);
    BOOST_TEST(g.getGlobalBin({0.25, 0}) == 1u);
    BOOST_TEST(g.getGlobalBin({0.5, 0}) == 2u);
    BOOST_TEST(g.getGlobalBin({0.75, 0}) == 3u);
    BOOST_TEST(g.getGlobalBin({1, 0}) == 4u);
    BOOST_TEST(g.getGlobalBin({0, 0.5}) == 5u);
    BOOST_TEST(g.getGlobalBin({0.25, 0.5}) == 6u);
    BOOST_TEST(g.getGlobalBin({0.5, 0.5}) == 7u);
    BOOST_TEST(g.getGlobalBin({0.75, 0.5}) == 8u);
    BOOST_TEST(g.getGlobalBin({1, 0.5}) == 9u);
    BOOST_TEST(g.getGlobalBin({0, 3}) == 10u);
    BOOST_TEST(g.getGlobalBin({0.25, 3}) == 11u);
    BOOST_TEST(g.getGlobalBin({0.5, 3}) == 12u);
    BOOST_TEST(g.getGlobalBin({0.75, 3}) == 13u);
    BOOST_TEST(g.getGlobalBin({1, 3}) == 14u);

    // test some arbitrary points
    BOOST_TEST(g.getGlobalBin({1.2, 0.3}) == 4u);
    BOOST_TEST(g.getGlobalBin({0.2, 1.3}) == 5u);
    BOOST_TEST(g.getGlobalBin({0.9, 1.8}) == 8u);
    BOOST_TEST(g.getGlobalBin({0.7, 2.1}) == 7u);
    BOOST_TEST(g.getGlobalBin({0.4, 0.3}) == 1u);

    // some outside points
    BOOST_TEST(g.getGlobalBin({-3, 2}) == 5u);
    BOOST_TEST(g.getGlobalBin({8, 1}) == 9u);
    BOOST_TEST(g.getGlobalBin({0.1, -3}) == 0u);
    BOOST_TEST(g.getGlobalBin({0.8, 11}) == 13u);
    BOOST_TEST(g.getGlobalBin({-2, -3}) == 0u);
    BOOST_TEST(g.getGlobalBin({-2, 7}) == 10u);
    BOOST_TEST(g.getGlobalBin({12, -1}) == 4u);
    BOOST_TEST(g.getGlobalBin({12, 11}) == 14u);
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
    BOOST_TEST(g.size() == 15u);

    // test some arbitrary points
    BOOST_TEST(g.at(Point({1.2, 0.3})) == 0.);
    BOOST_TEST(g.at(Point({2.2, 1.3})) == 1.);
    BOOST_TEST(g.at(Point({4.9, 1.8})) == 8.);
    BOOST_TEST(g.at(Point({3.7, 2.1})) == 7.);
    BOOST_TEST(g.at(Point({0.4, 2.3})) == 5.);
  }

}  // namespace Test

}  // namespace Acts
