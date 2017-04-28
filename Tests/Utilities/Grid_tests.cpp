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

  BOOST_AUTO_TEST_CASE(grid_test_1d_equidistant)
  {
    typedef std::array<double, 1> Point;
    typedef std::array<size_t, 1> indices;
    EquidistantAxis a(0.0, 4.0, 4u);
    Grid<double, EquidistantAxis> g(std::make_tuple(std::move(a)));

    // test general properties
    BOOST_TEST(g.size() == 6u);

    // global bin index
    BOOST_TEST(g.getGlobalBinIndex(Point({-0.3})) == 0u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-0.})) == 1u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.})) == 1u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.7})) == 1u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1})) == 2u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1.2})) == 2u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2.})) == 3u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2.7})) == 3u);
    BOOST_TEST(g.getGlobalBinIndex(Point({3.})) == 4u);
    BOOST_TEST(g.getGlobalBinIndex(Point({3.9999})) == 4u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4.})) == 5u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4.98})) == 5u);

    // global bin index -> local bin indices
    BOOST_TEST(g.getLocalBinIndices(0) == indices({0}));
    BOOST_TEST(g.getLocalBinIndices(1) == indices({1}));
    BOOST_TEST(g.getLocalBinIndices(2) == indices({2}));
    BOOST_TEST(g.getLocalBinIndices(3) == indices({3}));
    BOOST_TEST(g.getLocalBinIndices(4) == indices({4}));
    BOOST_TEST(g.getLocalBinIndices(5) == indices({5}));

    // local bin indices -> global bin index
    BOOST_TEST(g.getGlobalBinIndex({0}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({1}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({2}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({3}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({4}) == 4u);
    BOOST_TEST(g.getGlobalBinIndex({5}) == 5u);

    BOOST_TEST(g.getLocalBinIndices(g.getGlobalBinIndex(Point({2.7})))
               == indices({3}));

    // initialize grid
    for (size_t bin = 0; bin < g.size(); ++bin) g.at(bin) = bin;

    // test some bin centers
    BOOST_TEST(g.getBinCenter({1}) == Point({0.5}));
    BOOST_TEST(g.getBinCenter({2}) == Point({1.5}));
    BOOST_TEST(g.getBinCenter({3}) == Point({2.5}));
    BOOST_TEST(g.getBinCenter({4}) == Point({3.5}));

    // test some lower-left bin edges
    BOOST_TEST(g.getLowerLeftBinEdge({1}) == Point({0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({2}) == Point({1.}));
    BOOST_TEST(g.getLowerLeftBinEdge({3}) == Point({2.}));
    BOOST_TEST(g.getLowerLeftBinEdge({4}) == Point({3.}));

    // test some upper right-bin edges
    BOOST_TEST(g.getUpperRightBinEdge({1}) == Point({1.}));
    BOOST_TEST(g.getUpperRightBinEdge({2}) == Point({2.}));
    BOOST_TEST(g.getUpperRightBinEdge({3}) == Point({3.}));
    BOOST_TEST(g.getUpperRightBinEdge({4}) == Point({4.}));

    // consistency of access
    const auto& point     = Point({0.7});
    size_t      globalBin = g.getGlobalBinIndex(point);
    indices     localBins = g.getLocalBinIndices(globalBin);

    BOOST_TEST(g.at(point) == g.at(globalBin));
    BOOST_TEST(g.at(point) == g.at(localBins));
  }

  BOOST_AUTO_TEST_CASE(grid_test_2d_equidistant)
  {
    typedef std::array<double, 2> Point;
    typedef std::array<size_t, 2> indices;
    EquidistantAxis a(0.0, 4.0, 4u);
    EquidistantAxis b(0.0, 3.0, 3u);
    Grid<double, EquidistantAxis, EquidistantAxis> g(
        std::make_tuple(std::move(a), std::move(b)));

    // test general properties
    BOOST_TEST(g.size() == 30u);

    // global bin index
    BOOST_TEST(g.getGlobalBinIndex(Point({-1, -1})) == 0u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, -1})) == 1u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, -1})) == 2u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, -1})) == 3u);
    BOOST_TEST(g.getGlobalBinIndex(Point({3, -1})) == 4u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4, -1})) == 5u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-1, 0})) == 6u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0})) == 7u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0})) == 8u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 0})) == 9u);
    BOOST_TEST(g.getGlobalBinIndex(Point({3, 0})) == 10u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4, 0})) == 11u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-1, 1})) == 12u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 1})) == 13u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 1})) == 14u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 1})) == 15u);
    BOOST_TEST(g.getGlobalBinIndex(Point({3, 1})) == 16u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4, 1})) == 17u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-1, 2})) == 18u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 2})) == 19u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 2})) == 20u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 2})) == 21u);
    BOOST_TEST(g.getGlobalBinIndex(Point({3, 2})) == 22u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4, 2})) == 23u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-1, 3})) == 24u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 3})) == 25u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 3})) == 26u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 3})) == 27u);
    BOOST_TEST(g.getGlobalBinIndex(Point({3, 3})) == 28u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4, 3})) == 29u);

    // test some arbitrary points
    BOOST_TEST(g.getGlobalBinIndex(Point({1.2, 0.3})) == 8u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2.2, 3.3})) == 27u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.9, 1.8})) == 13u);
    BOOST_TEST(g.getGlobalBinIndex(Point({3.7, 3.1})) == 28u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1.4, 2.3})) == 20u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-3, 3})) == 24u);
    BOOST_TEST(g.getGlobalBinIndex(Point({8, 1})) == 17u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, -3})) == 2u);
    BOOST_TEST(g.getGlobalBinIndex(Point({3, 11})) == 28u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-2, -3})) == 0u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-2, 7})) == 24u);
    BOOST_TEST(g.getGlobalBinIndex(Point({12, -1})) == 5u);
    BOOST_TEST(g.getGlobalBinIndex(Point({12, 11})) == 29u);

    // global bin index -> local bin indices
    BOOST_TEST(g.getLocalBinIndices(0) == indices({0, 0}));
    BOOST_TEST(g.getLocalBinIndices(1) == indices({1, 0}));
    BOOST_TEST(g.getLocalBinIndices(2) == indices({2, 0}));
    BOOST_TEST(g.getLocalBinIndices(3) == indices({3, 0}));
    BOOST_TEST(g.getLocalBinIndices(4) == indices({4, 0}));
    BOOST_TEST(g.getLocalBinIndices(5) == indices({5, 0}));
    BOOST_TEST(g.getLocalBinIndices(6) == indices({0, 1}));
    BOOST_TEST(g.getLocalBinIndices(7) == indices({1, 1}));
    BOOST_TEST(g.getLocalBinIndices(8) == indices({2, 1}));
    BOOST_TEST(g.getLocalBinIndices(9) == indices({3, 1}));
    BOOST_TEST(g.getLocalBinIndices(10) == indices({4, 1}));
    BOOST_TEST(g.getLocalBinIndices(11) == indices({5, 1}));
    BOOST_TEST(g.getLocalBinIndices(12) == indices({0, 2}));
    BOOST_TEST(g.getLocalBinIndices(13) == indices({1, 2}));
    BOOST_TEST(g.getLocalBinIndices(14) == indices({2, 2}));
    BOOST_TEST(g.getLocalBinIndices(15) == indices({3, 2}));
    BOOST_TEST(g.getLocalBinIndices(16) == indices({4, 2}));
    BOOST_TEST(g.getLocalBinIndices(17) == indices({5, 2}));
    BOOST_TEST(g.getLocalBinIndices(18) == indices({0, 3}));
    BOOST_TEST(g.getLocalBinIndices(19) == indices({1, 3}));
    BOOST_TEST(g.getLocalBinIndices(20) == indices({2, 3}));
    BOOST_TEST(g.getLocalBinIndices(21) == indices({3, 3}));
    BOOST_TEST(g.getLocalBinIndices(22) == indices({4, 3}));
    BOOST_TEST(g.getLocalBinIndices(23) == indices({5, 3}));
    BOOST_TEST(g.getLocalBinIndices(24) == indices({0, 4}));
    BOOST_TEST(g.getLocalBinIndices(25) == indices({1, 4}));
    BOOST_TEST(g.getLocalBinIndices(26) == indices({2, 4}));
    BOOST_TEST(g.getLocalBinIndices(27) == indices({3, 4}));
    BOOST_TEST(g.getLocalBinIndices(28) == indices({4, 4}));
    BOOST_TEST(g.getLocalBinIndices(29) == indices({5, 4}));

    // local bin indices -> global bin index
    BOOST_TEST(g.getGlobalBinIndex({0, 0}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({3, 0}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({4, 0}) == 4u);
    BOOST_TEST(g.getGlobalBinIndex({5, 0}) == 5u);
    BOOST_TEST(g.getGlobalBinIndex({0, 1}) == 6u);
    BOOST_TEST(g.getGlobalBinIndex({1, 1}) == 7u);
    BOOST_TEST(g.getGlobalBinIndex({2, 1}) == 8u);
    BOOST_TEST(g.getGlobalBinIndex({3, 1}) == 9u);
    BOOST_TEST(g.getGlobalBinIndex({4, 1}) == 10u);
    BOOST_TEST(g.getGlobalBinIndex({5, 1}) == 11u);
    BOOST_TEST(g.getGlobalBinIndex({0, 2}) == 12u);
    BOOST_TEST(g.getGlobalBinIndex({1, 2}) == 13u);
    BOOST_TEST(g.getGlobalBinIndex({2, 2}) == 14u);
    BOOST_TEST(g.getGlobalBinIndex({3, 2}) == 15u);
    BOOST_TEST(g.getGlobalBinIndex({4, 2}) == 16u);
    BOOST_TEST(g.getGlobalBinIndex({5, 2}) == 17u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3}) == 18u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3}) == 19u);
    BOOST_TEST(g.getGlobalBinIndex({2, 3}) == 20u);
    BOOST_TEST(g.getGlobalBinIndex({3, 3}) == 21u);
    BOOST_TEST(g.getGlobalBinIndex({4, 3}) == 22u);
    BOOST_TEST(g.getGlobalBinIndex({5, 3}) == 23u);
    BOOST_TEST(g.getGlobalBinIndex({0, 4}) == 24u);
    BOOST_TEST(g.getGlobalBinIndex({1, 4}) == 25u);
    BOOST_TEST(g.getGlobalBinIndex({2, 4}) == 26u);
    BOOST_TEST(g.getGlobalBinIndex({3, 4}) == 27u);
    BOOST_TEST(g.getGlobalBinIndex({4, 4}) == 28u);
    BOOST_TEST(g.getGlobalBinIndex({5, 4}) == 29u);

    BOOST_TEST(g.getLocalBinIndices(g.getGlobalBinIndex(Point({1.2, 0.7})))
               == indices({2, 1}));

    // test some bin centers
    BOOST_TEST(g.getBinCenter({1, 1}) == Point({0.5, 0.5}));
    BOOST_TEST(g.getBinCenter({2, 3}) == Point({1.5, 2.5}));
    BOOST_TEST(g.getBinCenter({3, 1}) == Point({2.5, 0.5}));
    BOOST_TEST(g.getBinCenter({4, 2}) == Point({3.5, 1.5}));

    // test some lower-left bin edges
    BOOST_TEST(g.getLowerLeftBinEdge({1, 1}) == Point({0., 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({2, 3}) == Point({1., 2.}));
    BOOST_TEST(g.getLowerLeftBinEdge({3, 1}) == Point({2., 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({4, 2}) == Point({3., 1.}));

    // test some upper right-bin edges
    BOOST_TEST(g.getUpperRightBinEdge({1, 1}) == Point({1., 1.}));
    BOOST_TEST(g.getUpperRightBinEdge({2, 3}) == Point({2., 3.}));
    BOOST_TEST(g.getUpperRightBinEdge({3, 1}) == Point({3., 1.}));
    BOOST_TEST(g.getUpperRightBinEdge({4, 2}) == Point({4., 2.}));

    // initialize grid
    for (size_t bin = 0; bin < g.size(); ++bin) g.at(bin) = bin;

    // consistency of access
    const auto& point     = Point({0.7, 1.3});
    size_t      globalBin = g.getGlobalBinIndex(point);
    indices     localBins = g.getLocalBinIndices(globalBin);

    BOOST_TEST(g.at(point) == g.at(globalBin));
    BOOST_TEST(g.at(point) == g.at(localBins));
  }

  BOOST_AUTO_TEST_CASE(grid_test_3d_equidistant)
  {
    typedef std::array<double, 3> Point;
    typedef std::array<size_t, 3> indices;
    EquidistantAxis a(0.0, 2.0, 2u);
    EquidistantAxis b(0.0, 3.0, 3u);
    EquidistantAxis c(0.0, 2.0, 2u);
    Grid<double, EquidistantAxis, EquidistantAxis, EquidistantAxis> g(
        std::make_tuple(std::move(a), std::move(b), std::move(c)));

    // test general properties
    BOOST_TEST(g.size() == 80u);

    // test grid points
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0, 0})) == 25u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0, 0})) == 26u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 0, 0})) == 27u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 1, 0})) == 29u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 1, 0})) == 30u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 1, 0})) == 31u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 2, 0})) == 33u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 2, 0})) == 34u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 2, 0})) == 35u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 3, 0})) == 37u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 3, 0})) == 38u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 3, 0})) == 39u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0, 1})) == 45u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0, 1})) == 46u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 0, 1})) == 47u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 1, 1})) == 49u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 1, 1})) == 50u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 1, 1})) == 51u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 2, 1})) == 53u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 2, 1})) == 54u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 2, 1})) == 55u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 3, 1})) == 57u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 3, 1})) == 58u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 3, 1})) == 59u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0, 2})) == 65u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0, 2})) == 66u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 0, 2})) == 67u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 1, 2})) == 69u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 1, 2})) == 70u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 1, 2})) == 71u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 2, 2})) == 73u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 2, 2})) == 74u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 2, 2})) == 75u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 3, 2})) == 77u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 3, 2})) == 78u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2, 3, 2})) == 79u);

    // global bin index -> local bin indices
    BOOST_TEST(g.getLocalBinIndices(0) == indices({0, 0, 0}));
    BOOST_TEST(g.getLocalBinIndices(1) == indices({1, 0, 0}));
    BOOST_TEST(g.getLocalBinIndices(2) == indices({2, 0, 0}));
    BOOST_TEST(g.getLocalBinIndices(3) == indices({3, 0, 0}));
    BOOST_TEST(g.getLocalBinIndices(4) == indices({0, 1, 0}));
    BOOST_TEST(g.getLocalBinIndices(5) == indices({1, 1, 0}));
    BOOST_TEST(g.getLocalBinIndices(6) == indices({2, 1, 0}));
    BOOST_TEST(g.getLocalBinIndices(7) == indices({3, 1, 0}));
    BOOST_TEST(g.getLocalBinIndices(24) == indices({0, 1, 1}));
    BOOST_TEST(g.getLocalBinIndices(25) == indices({1, 1, 1}));
    BOOST_TEST(g.getLocalBinIndices(26) == indices({2, 1, 1}));
    BOOST_TEST(g.getLocalBinIndices(27) == indices({3, 1, 1}));
    BOOST_TEST(g.getLocalBinIndices(52) == indices({0, 3, 2}));
    BOOST_TEST(g.getLocalBinIndices(53) == indices({1, 3, 2}));
    BOOST_TEST(g.getLocalBinIndices(54) == indices({2, 3, 2}));
    BOOST_TEST(g.getLocalBinIndices(55) == indices({3, 3, 2}));
    BOOST_TEST(g.getLocalBinIndices(60) == indices({0, 0, 3}));
    BOOST_TEST(g.getLocalBinIndices(61) == indices({1, 0, 3}));
    BOOST_TEST(g.getLocalBinIndices(62) == indices({2, 0, 3}));
    BOOST_TEST(g.getLocalBinIndices(63) == indices({3, 0, 3}));
    BOOST_TEST(g.getLocalBinIndices(76) == indices({0, 4, 3}));
    BOOST_TEST(g.getLocalBinIndices(77) == indices({1, 4, 3}));
    BOOST_TEST(g.getLocalBinIndices(78) == indices({2, 4, 3}));
    BOOST_TEST(g.getLocalBinIndices(79) == indices({3, 4, 3}));

    // local bin indices -> global bin index
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 0}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 0}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0, 0}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({3, 0, 0}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({0, 1, 0}) == 4u);
    BOOST_TEST(g.getGlobalBinIndex({1, 1, 0}) == 5u);
    BOOST_TEST(g.getGlobalBinIndex({2, 1, 0}) == 6u);
    BOOST_TEST(g.getGlobalBinIndex({3, 1, 0}) == 7u);
    BOOST_TEST(g.getGlobalBinIndex({0, 1, 1}) == 24u);
    BOOST_TEST(g.getGlobalBinIndex({1, 1, 1}) == 25u);
    BOOST_TEST(g.getGlobalBinIndex({2, 1, 1}) == 26u);
    BOOST_TEST(g.getGlobalBinIndex({3, 1, 1}) == 27u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3, 2}) == 52u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3, 2}) == 53u);
    BOOST_TEST(g.getGlobalBinIndex({2, 3, 2}) == 54u);
    BOOST_TEST(g.getGlobalBinIndex({3, 3, 2}) == 55u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 3}) == 60u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 3}) == 61u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0, 3}) == 62u);
    BOOST_TEST(g.getGlobalBinIndex({3, 0, 3}) == 63u);
    BOOST_TEST(g.getGlobalBinIndex({0, 4, 3}) == 76u);
    BOOST_TEST(g.getGlobalBinIndex({1, 4, 3}) == 77u);
    BOOST_TEST(g.getGlobalBinIndex({2, 4, 3}) == 78u);
    BOOST_TEST(g.getGlobalBinIndex({3, 4, 3}) == 79u);

    BOOST_TEST(g.getLocalBinIndices(g.getGlobalBinIndex(Point({1.2, 0.7, 1.4})))
               == indices({2, 1, 2}));

    // test some bin centers
    BOOST_TEST(g.getBinCenter({1, 1, 1}) == Point({0.5, 0.5, 0.5}));
    BOOST_TEST(g.getBinCenter({2, 3, 2}) == Point({1.5, 2.5, 1.5}));
    BOOST_TEST(g.getBinCenter({1, 1, 2}) == Point({0.5, 0.5, 1.5}));
    BOOST_TEST(g.getBinCenter({2, 2, 1}) == Point({1.5, 1.5, 0.5}));

    // test some lower-left bin edges
    BOOST_TEST(g.getLowerLeftBinEdge({1, 1, 1}) == Point({0., 0., 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({2, 3, 2}) == Point({1., 2., 1.}));
    BOOST_TEST(g.getLowerLeftBinEdge({1, 1, 2}) == Point({0., 0., 1.}));
    BOOST_TEST(g.getLowerLeftBinEdge({2, 2, 1}) == Point({1., 1., 0.}));

    // test some upper right-bin edges
    BOOST_TEST(g.getUpperRightBinEdge({1, 1, 1}) == Point({1., 1., 1.}));
    BOOST_TEST(g.getUpperRightBinEdge({2, 3, 2}) == Point({2., 3., 2.}));
    BOOST_TEST(g.getUpperRightBinEdge({1, 1, 2}) == Point({1., 1., 2.}));
    BOOST_TEST(g.getUpperRightBinEdge({2, 2, 1}) == Point({2., 2., 1.}));

    // initialize grid
    for (size_t bin = 0; bin < g.size(); ++bin) g.at(bin) = bin;

    // consistency of access
    const auto& point     = Point({0.7, 2.3, 1.3});
    size_t      globalBin = g.getGlobalBinIndex(point);
    indices     localBins = g.getLocalBinIndices(globalBin);

    BOOST_TEST(g.at(point) == g.at(globalBin));
    BOOST_TEST(g.at(point) == g.at(localBins));
  }

  BOOST_AUTO_TEST_CASE(grid_test_1d_variable)
  {
    typedef std::array<double, 1> Point;
    typedef std::array<size_t, 1> indices;
    VariableAxis a({0.0, 1.0, 4.0});
    Grid<double, VariableAxis> g(std::make_tuple(std::move(a)));

    // test general properties
    BOOST_TEST(g.size() == 4u);

    // global bin index
    BOOST_TEST(g.getGlobalBinIndex(Point({-0.3})) == 0u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.})) == 1u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.7})) == 1u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1})) == 2u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1.2})) == 2u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2.7})) == 2u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4.})) == 3u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4.98})) == 3u);

    // global bin index -> local bin indices
    BOOST_TEST(g.getLocalBinIndices(0) == indices({0}));
    BOOST_TEST(g.getLocalBinIndices(1) == indices({1}));
    BOOST_TEST(g.getLocalBinIndices(2) == indices({2}));
    BOOST_TEST(g.getLocalBinIndices(3) == indices({3}));

    // local bin indices -> global bin index
    BOOST_TEST(g.getGlobalBinIndex({0}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({1}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({2}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({3}) == 3u);

    BOOST_TEST(g.getLocalBinIndices(g.getGlobalBinIndex(Point({0.8})))
               == indices({1}));

    // test some bin centers
    BOOST_TEST(g.getBinCenter({1}) == Point({0.5}));
    BOOST_TEST(g.getBinCenter({2}) == Point({2.5}));

    // test some lower-left bin edges
    BOOST_TEST(g.getLowerLeftBinEdge({1}) == Point({0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({2}) == Point({1.}));

    // test some upper right-bin edges
    BOOST_TEST(g.getUpperRightBinEdge({1}) == Point({1.}));
    BOOST_TEST(g.getUpperRightBinEdge({2}) == Point({4.}));

    // initialize grid
    for (size_t bin = 0; bin < g.size(); ++bin) g.at(bin) = bin;

    // consistency of access
    const auto& point     = Point({0.7});
    size_t      globalBin = g.getGlobalBinIndex(point);
    indices     localBins = g.getLocalBinIndices(globalBin);

    BOOST_TEST(g.at(point) == g.at(globalBin));
    BOOST_TEST(g.at(point) == g.at(localBins));
  }

  BOOST_AUTO_TEST_CASE(grid_test_2d_variable)
  {
    typedef std::array<double, 2> Point;
    typedef std::array<size_t, 2> indices;
    VariableAxis a({0.0, 1.0, 4.0});
    VariableAxis b({0.0, 0.5, 3.0});
    Grid<double, VariableAxis, VariableAxis> g(
        std::make_tuple(std::move(a), std::move(b)));

    // test general properties
    BOOST_TEST(g.size() == 16u);

    // test grid points
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0})) == 5u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0})) == 6u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4, 0})) == 7u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0.5})) == 9u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0.5})) == 10u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4, 0.5})) == 11u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 3})) == 13u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 3})) == 14u);
    BOOST_TEST(g.getGlobalBinIndex(Point({4, 3})) == 15u);

    // test some arbitrary points
    BOOST_TEST(g.getGlobalBinIndex(Point({1.2, 0.3})) == 6u);
    BOOST_TEST(g.getGlobalBinIndex(Point({2.2, 3.3})) == 14u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.9, 1.8})) == 9u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.7, 3.1})) == 13u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1.4, 2.3})) == 10u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-3, 2})) == 8u);
    BOOST_TEST(g.getGlobalBinIndex(Point({8, 1})) == 11u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, -3})) == 2u);
    BOOST_TEST(g.getGlobalBinIndex(Point({3, 11})) == 14u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-2, -3})) == 0u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-2, 7})) == 12u);
    BOOST_TEST(g.getGlobalBinIndex(Point({12, -1})) == 3u);
    BOOST_TEST(g.getGlobalBinIndex(Point({12, 11})) == 15u);

    // global bin index -> local bin indices
    BOOST_TEST(g.getLocalBinIndices(0) == indices({0, 0}));
    BOOST_TEST(g.getLocalBinIndices(1) == indices({1, 0}));
    BOOST_TEST(g.getLocalBinIndices(2) == indices({2, 0}));
    BOOST_TEST(g.getLocalBinIndices(3) == indices({3, 0}));
    BOOST_TEST(g.getLocalBinIndices(4) == indices({0, 1}));
    BOOST_TEST(g.getLocalBinIndices(5) == indices({1, 1}));
    BOOST_TEST(g.getLocalBinIndices(6) == indices({2, 1}));
    BOOST_TEST(g.getLocalBinIndices(7) == indices({3, 1}));
    BOOST_TEST(g.getLocalBinIndices(8) == indices({0, 2}));
    BOOST_TEST(g.getLocalBinIndices(9) == indices({1, 2}));
    BOOST_TEST(g.getLocalBinIndices(10) == indices({2, 2}));
    BOOST_TEST(g.getLocalBinIndices(11) == indices({3, 2}));
    BOOST_TEST(g.getLocalBinIndices(12) == indices({0, 3}));
    BOOST_TEST(g.getLocalBinIndices(13) == indices({1, 3}));
    BOOST_TEST(g.getLocalBinIndices(14) == indices({2, 3}));
    BOOST_TEST(g.getLocalBinIndices(15) == indices({3, 3}));

    // local bin indices -> global bin index
    BOOST_TEST(g.getGlobalBinIndex({0, 0}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({3, 0}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({0, 1}) == 4u);
    BOOST_TEST(g.getGlobalBinIndex({1, 1}) == 5u);
    BOOST_TEST(g.getGlobalBinIndex({2, 1}) == 6u);
    BOOST_TEST(g.getGlobalBinIndex({3, 1}) == 7u);
    BOOST_TEST(g.getGlobalBinIndex({0, 2}) == 8u);
    BOOST_TEST(g.getGlobalBinIndex({1, 2}) == 9u);
    BOOST_TEST(g.getGlobalBinIndex({2, 2}) == 10u);
    BOOST_TEST(g.getGlobalBinIndex({3, 2}) == 11u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3}) == 12u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3}) == 13u);
    BOOST_TEST(g.getGlobalBinIndex({2, 3}) == 14u);
    BOOST_TEST(g.getGlobalBinIndex({3, 3}) == 15u);

    BOOST_TEST(g.getLocalBinIndices(g.getGlobalBinIndex(Point({1.8, 3.2})))
               == indices({2, 3}));

    // test some bin centers
    BOOST_TEST(g.getBinCenter({1, 1}) == Point({0.5, 0.25}));
    BOOST_TEST(g.getBinCenter({1, 2}) == Point({0.5, 1.75}));
    BOOST_TEST(g.getBinCenter({2, 1}) == Point({2.5, 0.25}));
    BOOST_TEST(g.getBinCenter({2, 2}) == Point({2.5, 1.75}));

    // test some lower-left bin edges
    BOOST_TEST(g.getLowerLeftBinEdge({1, 1}) == Point({0., 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({1, 2}) == Point({0., 0.5}));
    BOOST_TEST(g.getLowerLeftBinEdge({2, 1}) == Point({1., 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({2, 2}) == Point({1., 0.5}));

    // test some upper right-bin edges
    BOOST_TEST(g.getUpperRightBinEdge({1, 1}) == Point({1., 0.5}));
    BOOST_TEST(g.getUpperRightBinEdge({1, 2}) == Point({1., 3.}));
    BOOST_TEST(g.getUpperRightBinEdge({2, 1}) == Point({4., 0.5}));
    BOOST_TEST(g.getUpperRightBinEdge({2, 2}) == Point({4., 3.}));

    // initialize grid
    for (size_t bin = 0; bin < g.size(); ++bin) g.at(bin) = bin;

    // consistency of access
    const auto& point     = Point({0.7, 1.3});
    size_t      globalBin = g.getGlobalBinIndex(point);
    indices     localBins = g.getLocalBinIndices(globalBin);

    BOOST_TEST(g.at(point) == g.at(globalBin));
    BOOST_TEST(g.at(point) == g.at(localBins));
  }

  BOOST_AUTO_TEST_CASE(grid_test_3d_variable)
  {
    typedef std::array<double, 3> Point;
    typedef std::array<size_t, 3> indices;
    VariableAxis a({0.0, 1.0});
    VariableAxis b({0.0, 0.5, 3.0});
    VariableAxis c({0.0, 0.5, 3.0, 3.3});
    Grid<double, VariableAxis, VariableAxis, VariableAxis> g(
        std::make_tuple(std::move(a), std::move(b), std::move(c)));

    // test general properties
    BOOST_TEST(g.size() == 60u);

    // test grid points
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0, 0})) == 16u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0, 0})) == 17u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0.5, 0})) == 19u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0.5, 0})) == 20u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 3, 0})) == 22u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 3, 0})) == 23u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0, 0.5})) == 28u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0, 0.5})) == 29u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0.5, 0.5})) == 31u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0.5, 0.5})) == 32u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 3, 0.5})) == 34u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 3, 0.5})) == 35u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0, 3})) == 40u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0, 3})) == 41u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0.5, 3})) == 43u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0.5, 3})) == 44u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 3, 3})) == 46u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 3, 3})) == 47u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0, 3.3})) == 52u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0, 3.3})) == 53u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0.5, 3.3})) == 55u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0.5, 3.3})) == 56u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 3, 3.3})) == 58u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 3, 3.3})) == 59u);

    // global bin index -> local bin indices
    BOOST_TEST(g.getLocalBinIndices(0) == indices({0, 0, 0}));
    BOOST_TEST(g.getLocalBinIndices(1) == indices({1, 0, 0}));
    BOOST_TEST(g.getLocalBinIndices(2) == indices({2, 0, 0}));
    BOOST_TEST(g.getLocalBinIndices(3) == indices({0, 1, 0}));
    BOOST_TEST(g.getLocalBinIndices(4) == indices({1, 1, 0}));
    BOOST_TEST(g.getLocalBinIndices(5) == indices({2, 1, 0}));
    BOOST_TEST(g.getLocalBinIndices(21) == indices({0, 3, 1}));
    BOOST_TEST(g.getLocalBinIndices(22) == indices({1, 3, 1}));
    BOOST_TEST(g.getLocalBinIndices(23) == indices({2, 3, 1}));
    BOOST_TEST(g.getLocalBinIndices(24) == indices({0, 0, 2}));
    BOOST_TEST(g.getLocalBinIndices(25) == indices({1, 0, 2}));
    BOOST_TEST(g.getLocalBinIndices(26) == indices({2, 0, 2}));
    BOOST_TEST(g.getLocalBinIndices(57) == indices({0, 3, 4}));
    BOOST_TEST(g.getLocalBinIndices(58) == indices({1, 3, 4}));
    BOOST_TEST(g.getLocalBinIndices(59) == indices({2, 3, 4}));

    // local bin indices -> global bin index
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 0}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 0}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0, 0}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({0, 1, 0}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({1, 1, 0}) == 4u);
    BOOST_TEST(g.getGlobalBinIndex({2, 1, 0}) == 5u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3, 1}) == 21u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3, 1}) == 22u);
    BOOST_TEST(g.getGlobalBinIndex({2, 3, 1}) == 23u);
    BOOST_TEST(g.getGlobalBinIndex({0, 0, 2}) == 24u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0, 2}) == 25u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0, 2}) == 26u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3, 4}) == 57u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3, 4}) == 58u);
    BOOST_TEST(g.getGlobalBinIndex({2, 3, 4}) == 59u);

    BOOST_TEST(g.getLocalBinIndices(g.getGlobalBinIndex(Point({1.8, 0.7, 3.2})))
               == indices({2, 2, 3}));

    // test some bin centers
    BOOST_TEST(g.getBinCenter({1, 1, 1}) == Point({0.5, 0.25, 0.25}));
    BOOST_TEST(g.getBinCenter({1, 1, 2}) == Point({0.5, 0.25, 1.75}));
    BOOST_TEST(g.getBinCenter({1, 1, 3}) == Point({0.5, 0.25, 3.15}));
    BOOST_TEST(g.getBinCenter({1, 2, 1}) == Point({0.5, 1.75, 0.25}));
    BOOST_TEST(g.getBinCenter({1, 2, 2}) == Point({0.5, 1.75, 1.75}));
    BOOST_TEST(g.getBinCenter({1, 2, 3}) == Point({0.5, 1.75, 3.15}));

    // test some lower-left bin edges
    BOOST_TEST(g.getLowerLeftBinEdge({1, 1, 1}) == Point({0., 0., 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({1, 1, 2}) == Point({0., 0., 0.5}));
    BOOST_TEST(g.getLowerLeftBinEdge({1, 1, 3}) == Point({0., 0., 3.}));
    BOOST_TEST(g.getLowerLeftBinEdge({1, 2, 1}) == Point({0., 0.5, 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({1, 2, 2}) == Point({0., 0.5, 0.5}));
    BOOST_TEST(g.getLowerLeftBinEdge({1, 2, 3}) == Point({0., 0.5, 3.}));

    // test some upper right-bin edges
    BOOST_TEST(g.getUpperRightBinEdge({1, 1, 1}) == Point({1., 0.5, 0.5}));
    BOOST_TEST(g.getUpperRightBinEdge({1, 1, 2}) == Point({1., 0.5, 3.}));
    BOOST_TEST(g.getUpperRightBinEdge({1, 1, 3}) == Point({1., 0.5, 3.3}));
    BOOST_TEST(g.getUpperRightBinEdge({1, 2, 1}) == Point({1., 3., 0.5}));
    BOOST_TEST(g.getUpperRightBinEdge({1, 2, 2}) == Point({1., 3., 3.}));
    BOOST_TEST(g.getUpperRightBinEdge({1, 2, 3}) == Point({1., 3., 3.3}));

    // initialize grid
    for (size_t bin = 0; bin < g.size(); ++bin) g.at(bin) = bin;

    // consistency of access
    const auto& point     = Point({0.7, 1.3, 3.7});
    size_t      globalBin = g.getGlobalBinIndex(point);
    indices     localBins = g.getLocalBinIndices(globalBin);

    BOOST_TEST(g.at(point) == g.at(globalBin));
    BOOST_TEST(g.at(point) == g.at(localBins));
  }

  BOOST_AUTO_TEST_CASE(grid_test_2d_mixed)
  {
    typedef std::array<double, 2> Point;
    typedef std::array<size_t, 2> indices;
    EquidistantAxis a(0.0, 1.0, 4u);
    VariableAxis    b({0.0, 0.5, 3.0});
    Grid<double, EquidistantAxis, VariableAxis> g(
        std::make_tuple(std::move(a), std::move(b)));

    // test general properties
    BOOST_TEST(g.size() == 24u);

    // test grid points
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0})) == 7u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.25, 0})) == 8u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.5, 0})) == 9u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.75, 0})) == 10u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0})) == 11u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 0.5})) == 13u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.25, 0.5})) == 14u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.5, 0.5})) == 15u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.75, 0.5})) == 16u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 0.5})) == 17u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0, 3})) == 19u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.25, 3})) == 20u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.5, 3})) == 21u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.75, 3})) == 22u);
    BOOST_TEST(g.getGlobalBinIndex(Point({1, 3})) == 23u);

    // test some arbitrary points
    BOOST_TEST(g.getGlobalBinIndex(Point({1.2, 0.3})) == 11u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.2, 1.3})) == 13u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.9, 1.8})) == 16u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.7, 2.1})) == 15u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.4, 0.3})) == 8u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-3, 2})) == 12u);
    BOOST_TEST(g.getGlobalBinIndex(Point({8, 1})) == 17u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.1, -3})) == 1u);
    BOOST_TEST(g.getGlobalBinIndex(Point({0.8, 11})) == 22u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-2, -3})) == 0u);
    BOOST_TEST(g.getGlobalBinIndex(Point({-2, 7})) == 18u);
    BOOST_TEST(g.getGlobalBinIndex(Point({12, -1})) == 5u);
    BOOST_TEST(g.getGlobalBinIndex(Point({12, 11})) == 23u);

    // global bin index -> local bin indices
    typedef std::array<size_t, 2> indices;
    BOOST_TEST(g.getLocalBinIndices(0) == indices({0, 0}));
    BOOST_TEST(g.getLocalBinIndices(1) == indices({1, 0}));
    BOOST_TEST(g.getLocalBinIndices(2) == indices({2, 0}));
    BOOST_TEST(g.getLocalBinIndices(3) == indices({3, 0}));
    BOOST_TEST(g.getLocalBinIndices(4) == indices({4, 0}));
    BOOST_TEST(g.getLocalBinIndices(5) == indices({5, 0}));
    BOOST_TEST(g.getLocalBinIndices(6) == indices({0, 1}));
    BOOST_TEST(g.getLocalBinIndices(7) == indices({1, 1}));
    BOOST_TEST(g.getLocalBinIndices(8) == indices({2, 1}));
    BOOST_TEST(g.getLocalBinIndices(9) == indices({3, 1}));
    BOOST_TEST(g.getLocalBinIndices(10) == indices({4, 1}));
    BOOST_TEST(g.getLocalBinIndices(11) == indices({5, 1}));
    BOOST_TEST(g.getLocalBinIndices(12) == indices({0, 2}));
    BOOST_TEST(g.getLocalBinIndices(13) == indices({1, 2}));
    BOOST_TEST(g.getLocalBinIndices(14) == indices({2, 2}));
    BOOST_TEST(g.getLocalBinIndices(15) == indices({3, 2}));
    BOOST_TEST(g.getLocalBinIndices(16) == indices({4, 2}));
    BOOST_TEST(g.getLocalBinIndices(17) == indices({5, 2}));
    BOOST_TEST(g.getLocalBinIndices(18) == indices({0, 3}));
    BOOST_TEST(g.getLocalBinIndices(19) == indices({1, 3}));
    BOOST_TEST(g.getLocalBinIndices(20) == indices({2, 3}));
    BOOST_TEST(g.getLocalBinIndices(21) == indices({3, 3}));
    BOOST_TEST(g.getLocalBinIndices(22) == indices({4, 3}));
    BOOST_TEST(g.getLocalBinIndices(23) == indices({5, 3}));

    // local bin indices -> global bin index
    BOOST_TEST(g.getGlobalBinIndex({0, 0}) == 0u);
    BOOST_TEST(g.getGlobalBinIndex({1, 0}) == 1u);
    BOOST_TEST(g.getGlobalBinIndex({2, 0}) == 2u);
    BOOST_TEST(g.getGlobalBinIndex({3, 0}) == 3u);
    BOOST_TEST(g.getGlobalBinIndex({4, 0}) == 4u);
    BOOST_TEST(g.getGlobalBinIndex({5, 0}) == 5u);
    BOOST_TEST(g.getGlobalBinIndex({0, 1}) == 6u);
    BOOST_TEST(g.getGlobalBinIndex({1, 1}) == 7u);
    BOOST_TEST(g.getGlobalBinIndex({2, 1}) == 8u);
    BOOST_TEST(g.getGlobalBinIndex({3, 1}) == 9u);
    BOOST_TEST(g.getGlobalBinIndex({4, 1}) == 10u);
    BOOST_TEST(g.getGlobalBinIndex({5, 1}) == 11u);
    BOOST_TEST(g.getGlobalBinIndex({0, 2}) == 12u);
    BOOST_TEST(g.getGlobalBinIndex({1, 2}) == 13u);
    BOOST_TEST(g.getGlobalBinIndex({2, 2}) == 14u);
    BOOST_TEST(g.getGlobalBinIndex({3, 2}) == 15u);
    BOOST_TEST(g.getGlobalBinIndex({4, 2}) == 16u);
    BOOST_TEST(g.getGlobalBinIndex({5, 2}) == 17u);
    BOOST_TEST(g.getGlobalBinIndex({0, 3}) == 18u);
    BOOST_TEST(g.getGlobalBinIndex({1, 3}) == 19u);
    BOOST_TEST(g.getGlobalBinIndex({2, 3}) == 20u);
    BOOST_TEST(g.getGlobalBinIndex({3, 3}) == 21u);
    BOOST_TEST(g.getGlobalBinIndex({4, 3}) == 22u);
    BOOST_TEST(g.getGlobalBinIndex({5, 3}) == 23u);

    BOOST_TEST(g.getLocalBinIndices(g.getGlobalBinIndex(Point({1.1, 1.7})))
               == indices({5, 2}));

    // test some bin centers
    BOOST_TEST(g.getBinCenter({1, 1}) == Point({0.125, 0.25}));
    BOOST_TEST(g.getBinCenter({1, 2}) == Point({0.125, 1.75}));
    BOOST_TEST(g.getBinCenter({2, 1}) == Point({0.375, 0.25}));
    BOOST_TEST(g.getBinCenter({2, 2}) == Point({0.375, 1.75}));
    BOOST_TEST(g.getBinCenter({3, 1}) == Point({0.625, 0.25}));
    BOOST_TEST(g.getBinCenter({3, 2}) == Point({0.625, 1.75}));
    BOOST_TEST(g.getBinCenter({4, 1}) == Point({0.875, 0.25}));
    BOOST_TEST(g.getBinCenter({4, 2}) == Point({0.875, 1.75}));

    // test some lower-left bin edges
    BOOST_TEST(g.getLowerLeftBinEdge({1, 1}) == Point({0., 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({1, 2}) == Point({0., 0.5}));
    BOOST_TEST(g.getLowerLeftBinEdge({2, 1}) == Point({0.25, 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({2, 2}) == Point({0.25, 0.5}));
    BOOST_TEST(g.getLowerLeftBinEdge({3, 1}) == Point({0.5, 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({3, 2}) == Point({0.5, 0.5}));
    BOOST_TEST(g.getLowerLeftBinEdge({4, 1}) == Point({0.75, 0.}));
    BOOST_TEST(g.getLowerLeftBinEdge({4, 2}) == Point({0.75, 0.5}));

    // test some upper-right bin edges
    BOOST_TEST(g.getUpperRightBinEdge({1, 1}) == Point({0.25, 0.5}));
    BOOST_TEST(g.getUpperRightBinEdge({1, 2}) == Point({0.25, 3.}));
    BOOST_TEST(g.getUpperRightBinEdge({2, 1}) == Point({0.5, 0.5}));
    BOOST_TEST(g.getUpperRightBinEdge({2, 2}) == Point({0.5, 3.}));
    BOOST_TEST(g.getUpperRightBinEdge({3, 1}) == Point({0.75, 0.5}));
    BOOST_TEST(g.getUpperRightBinEdge({3, 2}) == Point({0.75, 3.}));
    BOOST_TEST(g.getUpperRightBinEdge({4, 1}) == Point({1., 0.5}));
    BOOST_TEST(g.getUpperRightBinEdge({4, 2}) == Point({1., 3.}));

    // initialize grid
    for (size_t bin = 0; bin < g.size(); ++bin) g.at(bin) = bin;

    // consistency of access
    const auto& point     = Point({1.3, 3.7});
    size_t      globalBin = g.getGlobalBinIndex(point);
    indices     localBins = g.getLocalBinIndices(globalBin);

    BOOST_TEST(g.at(point) == g.at(globalBin));
    BOOST_TEST(g.at(point) == g.at(localBins));
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

  BOOST_AUTO_TEST_CASE(grid_interpolation)
  {
    typedef std::array<double, 3> Point;
    EquidistantAxis a(1.0, 3.0, 2u);
    EquidistantAxis b(1.0, 5.0, 2u);
    EquidistantAxis c(1.0, 7.0, 2u);
    Grid<double, EquidistantAxis, EquidistantAxis, EquidistantAxis> g(
        std::make_tuple(std::move(a), std::move(b), std::move(c)));

    g.at(Point({1., 1., 1.})) = 10.;
    g.at(Point({2., 1., 1.})) = 20.;
    g.at(Point({1., 3., 1.})) = 30.;
    g.at(Point({2., 3., 1.})) = 40.;
    g.at(Point({1., 1., 4.})) = 50.;
    g.at(Point({2., 1., 4.})) = 60.;
    g.at(Point({1., 3., 4.})) = 70.;
    g.at(Point({2., 3., 4.})) = 80.;

    BOOST_TEST(g.interpolate(Point({1., 1., 1.})) == 10.);
    BOOST_TEST(g.interpolate(Point({2., 1., 1.})) == 20.);
    BOOST_TEST(g.interpolate(Point({1., 3., 1.})) == 30.);
    BOOST_TEST(g.interpolate(Point({2., 3., 1.})) == 40.);
    BOOST_TEST(g.interpolate(Point({1., 1., 4.})) == 50.);
    BOOST_TEST(g.interpolate(Point({2., 1., 4.})) == 60.);
    BOOST_TEST(g.interpolate(Point({1., 3., 4.})) == 70.);
    BOOST_TEST(g.interpolate(Point({2., 3., 4.})) == 80.);
    BOOST_TEST(g.interpolate(Point({1.5, 1., 1.})) == 15.);
    BOOST_TEST(g.interpolate(Point({1.5, 3., 1.})) == 35.);
    BOOST_TEST(g.interpolate(Point({1., 2., 1.})) == 20.);
    BOOST_TEST(g.interpolate(Point({2., 2., 1.})) == 30.);
    BOOST_TEST(g.interpolate(Point({1.5, 1., 4.})) == 55.);
    BOOST_TEST(g.interpolate(Point({1.5, 3., 4.})) == 75.);
    BOOST_TEST(g.interpolate(Point({1., 2., 4.})) == 60.);
    BOOST_TEST(g.interpolate(Point({2., 2., 4.})) == 70.);
    BOOST_TEST(g.interpolate(Point({1., 1., 2.5})) == 30.);
    BOOST_TEST(g.interpolate(Point({1., 3., 2.5})) == 50.);
    BOOST_TEST(g.interpolate(Point({2., 1., 2.5})) == 40.);
    BOOST_TEST(g.interpolate(Point({2., 3., 2.5})) == 60.);
    BOOST_TEST(g.interpolate(Point({1.5, 2., 2.5})) == 360. / 8);
    BOOST_TEST(g.interpolate(Point({1.3, 2.1, 1.6})) == 32.);
  }
}  // namespace Test

}  // namespace Acts
