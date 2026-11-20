// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/MultiAxis.hpp"
#include "ActsTests/CommonHelpers/Assertions.hpp"
#include "ActsTests/CommonHelpers/FloatComparisons.hpp"

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

using namespace Acts;
using namespace Acts::detail;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(MultiAxisTests)

BOOST_AUTO_TEST_CASE(test_1d_equidistant) {
  using Point = std::array<double, 1>;
  using MultiIndex = std::array<std::size_t, 1>;

  const Axis a(0.0, 4.0, 4u);

  const MultiAxis ma(a);

  // test general properties
  BOOST_CHECK_EQUAL(ma.getNAxes(), 1u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(0), 4u);
  BOOST_CHECK_EQUAL(ma.getNTotalBins(), 4u);

  // flat bin index
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-0.3}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-0.}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.7}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1.2}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2.}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2.7}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3.}), 4u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3.9999}), 4u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({4.}), 5u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({4.98}), 5u);

  // flat bin index -> multi bin indices
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(0) == MultiIndex{0});
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(1) == MultiIndex{1});
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(2) == MultiIndex{2});
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(3) == MultiIndex{3});
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(4) == MultiIndex{4});
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(5) == MultiIndex{5});

  // multi bin indices -> flat bin index
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({4}), 4u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({5}), 5u);

  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(ma.getGlobalBinFromPoint({2.7})) ==
              MultiIndex{3});

  // inside checks
  BOOST_CHECK(!ma.isInside({-2.}));
  BOOST_CHECK(ma.isInside({0.}));
  BOOST_CHECK(ma.isInside({2.5}));
  BOOST_CHECK(!ma.isInside({4.}));
  BOOST_CHECK(!ma.isInside({6.}));

  // test some bin centers
  CHECK_CLOSE_ABS(ma.getBinCenter({1}), Point{0.5}, 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({2}), Point{1.5}, 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({3}), Point{2.5}, 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({4}), Point{3.5}, 1e-6);

  // test some lower-left bin edges
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1}), Point{0.}, 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({2}), Point{1.}, 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({3}), Point{2.}, 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({4}), Point{3.}, 1e-6);

  // test some upper right-bin edges
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1}), Point{1.}, 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({2}), Point{2.}, 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({3}), Point{3.}, 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({4}), Point{4.}, 1e-6);
}

BOOST_AUTO_TEST_CASE(test_2d_equidistant) {
  using Point = std::array<double, 2>;
  using MultiIndex = std::array<std::size_t, 2>;

  const Axis a(0.0, 4.0, 4u);
  const Axis b(0.0, 3.0, 3u);

  const MultiAxis ma(a, b);

  // test general properties
  BOOST_CHECK_EQUAL(ma.getNAxes(), 2u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(0), 4u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(1), 3u);
  BOOST_CHECK_EQUAL(ma.getNTotalBins(), 12u);

  // flat bin index
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-1, -1}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-1, 0}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-1, 1}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-1, 2}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-1, 3}), 4u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, -1}), 5u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0}), 6u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 1}), 7u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 2}), 8u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 3}), 9u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, -1}), 10u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0}), 11u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 1}), 12u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 2}), 13u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 3}), 14u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, -1}), 15u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 0}), 16u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 1}), 17u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 2}), 18u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 3}), 19u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3, -1}), 20u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3, 0}), 21u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3, 1}), 22u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3, 2}), 23u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3, 3}), 24u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({4, -1}), 25u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({4, 0}), 26u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({4, 1}), 27u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({4, 2}), 28u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({4, 3}), 29u);

  // test some arbitrary points
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1.2, 0.3}), 11u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2.2, 3.3}), 19u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.9, 1.8}), 7u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3.7, 3.1}), 24u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1.4, 2.3}), 13u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-3, 3}), 4u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({8, 1}), 27u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, -3}), 10u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3, 11}), 24u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-2, -3}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-2, 7}), 04u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({12, -1}), 25u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({12, 11}), 29u);

  // flat bin index -> multi bin indices
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(0) == (MultiIndex{0, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(1) == (MultiIndex{0, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(2) == (MultiIndex{0, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(3) == (MultiIndex{0, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(4) == (MultiIndex{0, 4}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(5) == (MultiIndex{1, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(6) == (MultiIndex{1, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(7) == (MultiIndex{1, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(8) == (MultiIndex{1, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(9) == (MultiIndex{1, 4}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(10) == (MultiIndex{2, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(11) == (MultiIndex{2, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(12) == (MultiIndex{2, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(13) == (MultiIndex{2, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(14) == (MultiIndex{2, 4}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(15) == (MultiIndex{3, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(16) == (MultiIndex{3, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(17) == (MultiIndex{3, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(18) == (MultiIndex{3, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(19) == (MultiIndex{3, 4}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(20) == (MultiIndex{4, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(21) == (MultiIndex{4, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(22) == (MultiIndex{4, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(23) == (MultiIndex{4, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(24) == (MultiIndex{4, 4}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(25) == (MultiIndex{5, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(26) == (MultiIndex{5, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(27) == (MultiIndex{5, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(28) == (MultiIndex{5, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(29) == (MultiIndex{5, 4}));

  // local bin indices -> global bin index
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 0}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 1}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 2}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 3}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 4}), 4u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 0}), 5u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 1}), 6u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 2}), 7u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 3}), 8u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 4}), 9u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 0}), 10u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 1}), 11u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 2}), 12u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 3}), 13u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 4}), 14u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 0}), 15u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 1}), 16u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 2}), 17u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 3}), 18u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 4}), 19u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({4, 0}), 20u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({4, 1}), 21u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({4, 2}), 22u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({4, 3}), 23u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({4, 4}), 24u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({5, 0}), 25u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({5, 1}), 26u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({5, 2}), 27u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({5, 3}), 28u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({5, 4}), 29u);

  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(
                  ma.getGlobalBinFromPoint({1.2, 0.7})) == (MultiIndex{2, 1}));

  // inside checks
  BOOST_CHECK(!ma.isInside({-2., -1}));
  BOOST_CHECK(!ma.isInside({-2., 1.}));
  BOOST_CHECK(!ma.isInside({-2., 5.}));
  BOOST_CHECK(!ma.isInside({1., -1.}));
  BOOST_CHECK(!ma.isInside({6., -1.}));
  BOOST_CHECK(ma.isInside({0.5, 1.3}));
  BOOST_CHECK(!ma.isInside({4., -1.}));
  BOOST_CHECK(!ma.isInside({4., 0.3}));
  BOOST_CHECK(!ma.isInside({4., 3.}));
  BOOST_CHECK(!ma.isInside({-1., 3.}));
  BOOST_CHECK(!ma.isInside({2., 3.}));
  BOOST_CHECK(!ma.isInside({5., 3.}));

  // test some bin centers
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 1}), (Point{0.5, 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({2, 3}), (Point{1.5, 2.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({3, 1}), (Point{2.5, 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({4, 2}), (Point{3.5, 1.5}), 1e-6);

  // test some lower-left bin edges
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 1}), (Point{0., 0.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({2, 3}), (Point{1., 2.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({3, 1}), (Point{2., 0.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({4, 2}), (Point{3., 1.}), 1e-6);

  // test some upper right-bin edges
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1, 1}), (Point{1., 1.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({2, 3}), (Point{2., 3.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({3, 1}), (Point{3., 1.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({4, 2}), (Point{4., 2.}), 1e-6);
}

BOOST_AUTO_TEST_CASE(test_3d_equidistant) {
  using Point = std::array<double, 3>;
  using MultiIndex = std::array<std::size_t, 3>;

  const Axis a(0.0, 2.0, 2u);
  const Axis b(0.0, 3.0, 3u);
  const Axis c(0.0, 2.0, 2u);

  const MultiAxis ma(a, b, c);

  // test general properties
  BOOST_CHECK_EQUAL(ma.getNAxes(), 3u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(0), 2u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(1), 3u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(2), 2u);
  BOOST_CHECK_EQUAL(ma.getNTotalBins(), 12u);

  // test grid points
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0, 0}), 25u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0, 1}), 26u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0, 2}), 27u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 1, 0}), 29u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 1, 1}), 30u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 1, 2}), 31u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 2, 0}), 33u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 2, 1}), 34u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 2, 2}), 35u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 3, 0}), 37u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 3, 1}), 38u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 3, 2}), 39u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0, 0}), 45u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0, 1}), 46u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0, 2}), 47u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 1, 0}), 49u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 1, 1}), 50u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 2, 0}), 53u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 2, 1}), 54u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 2, 2}), 55u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 3, 0}), 57u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 3, 1}), 58u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 3, 2}), 59u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 0, 0}), 65u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 0, 1}), 66u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 0, 2}), 67u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 1, 0}), 69u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 1, 1}), 70u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 1, 2}), 71u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 2, 0}), 73u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 2, 1}), 74u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 2, 2}), 75u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 3, 0}), 77u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 3, 1}), 78u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, 3, 2}), 79u);

  // flat bin index -> multi bin indices
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(0) == (MultiIndex{0, 0, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(1) == (MultiIndex{0, 0, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(2) == (MultiIndex{0, 0, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(3) == (MultiIndex{0, 0, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(4) == (MultiIndex{0, 1, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(5) == (MultiIndex{0, 1, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(6) == (MultiIndex{0, 1, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(7) == (MultiIndex{0, 1, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(24) == (MultiIndex{1, 1, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(25) == (MultiIndex{1, 1, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(26) == (MultiIndex{1, 1, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(27) == (MultiIndex{1, 1, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(52) == (MultiIndex{2, 3, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(53) == (MultiIndex{2, 3, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(54) == (MultiIndex{2, 3, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(55) == (MultiIndex{2, 3, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(60) == (MultiIndex{3, 0, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(61) == (MultiIndex{3, 0, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(62) == (MultiIndex{3, 0, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(63) == (MultiIndex{3, 0, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(76) == (MultiIndex{3, 4, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(77) == (MultiIndex{3, 4, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(78) == (MultiIndex{3, 4, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(79) == (MultiIndex{3, 4, 3}));

  // multi bin indices -> flat bin index
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 0, 0}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 0, 1}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 0, 2}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 0, 3}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 1, 0}), 4u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 1, 1}), 5u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 1, 2}), 6u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 1, 3}), 7u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 1, 0}), 24u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 1, 1}), 25u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 1, 2}), 26u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 1, 3}), 27u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 3, 0}), 52u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 3, 1}), 53u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 3, 2}), 54u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 3, 3}), 55u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 0, 0}), 60u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 0, 1}), 61u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 0, 2}), 62u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 0, 3}), 63u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 4, 0}), 76u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 4, 1}), 77u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 4, 2}), 78u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 4, 3}), 79u);

  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(ma.getGlobalBinFromPoint(
                  {1.2, 0.7, 1.4})) == (MultiIndex{2, 1, 2}));

  // inside checks
  BOOST_CHECK(!ma.isInside({-2., -1, -2}));
  BOOST_CHECK(!ma.isInside({-2., 1., 0.}));
  BOOST_CHECK(!ma.isInside({-2., 5., -1}));
  BOOST_CHECK(!ma.isInside({1., -1., 1.}));
  BOOST_CHECK(!ma.isInside({6., -1., 4.}));
  BOOST_CHECK(ma.isInside({0.5, 1.3, 1.7}));
  BOOST_CHECK(!ma.isInside({2., -1., -0.4}));
  BOOST_CHECK(!ma.isInside({2., 0.3, 3.4}));
  BOOST_CHECK(!ma.isInside({2., 3., 0.8}));
  BOOST_CHECK(!ma.isInside({-1., 3., 5.}));
  BOOST_CHECK(!ma.isInside({2., 3., -1.}));
  BOOST_CHECK(!ma.isInside({5., 3., 0.5}));

  // test some bin centers
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 1, 1}), (Point{0.5, 0.5, 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({2, 3, 2}), (Point{1.5, 2.5, 1.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 1, 2}), (Point{0.5, 0.5, 1.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({2, 2, 1}), (Point{1.5, 1.5, 0.5}), 1e-6);

  // test some lower-left bin edges
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 1, 1}), (Point{0., 0., 0.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({2, 3, 2}), (Point{1., 2., 1.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 1, 2}), (Point{0., 0., 1.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({2, 2, 1}), (Point{1., 1., 0.}), 1e-6);

  // test some upper right-bin edges
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({{1, 1, 1}}), (Point{{1., 1., 1.}}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({{2, 3, 2}}), (Point{{2., 3., 2.}}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({{1, 1, 2}}), (Point{{1., 1., 2.}}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({{2, 2, 1}}), (Point{{2., 2., 1.}}),
                  1e-6);
}

BOOST_AUTO_TEST_CASE(test_1d_variable) {
  using Point = std::array<double, 1>;
  using MultiIndex = std::array<std::size_t, 1>;

  const Axis a({0.0, 1.0, 4.0});

  const MultiAxis ma(a);

  // test general properties
  BOOST_CHECK_EQUAL(ma.getNAxes(), 1u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(0), 2u);
  BOOST_CHECK_EQUAL(ma.getNTotalBins(), 2u);

  // flat bin index
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-0.3}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.7}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1.2}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2.7}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({4.}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({4.98}), 3u);

  // flat bin index -> multi bin indices
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(0) == MultiIndex{0});
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(1) == MultiIndex{1});
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(2) == MultiIndex{2});
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(3) == MultiIndex{3});

  // multi bin indices -> flat bin index
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3}), 3u);

  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(ma.getGlobalBinFromPoint({0.8})) ==
              MultiIndex{1});

  // inside checks
  BOOST_CHECK(!ma.isInside({-2.}));
  BOOST_CHECK(ma.isInside({0.}));
  BOOST_CHECK(ma.isInside({2.5}));
  BOOST_CHECK(!ma.isInside({4.}));
  BOOST_CHECK(!ma.isInside({6.}));

  // test some bin centers
  CHECK_CLOSE_ABS(ma.getBinCenter({1}), Point{0.5}, 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({2}), Point{2.5}, 1e-6);

  // test some lower-left bin edges
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1}), Point{0.}, 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({2}), Point{1.}, 1e-6);

  // test some upper right-bin edges
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1}), Point{1.}, 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({2}), Point{4.}, 1e-6);
}

BOOST_AUTO_TEST_CASE(test_2d_variable) {
  using Point = std::array<double, 2>;
  using MultiIndex = std::array<std::size_t, 2>;

  const Axis a({0.0, 0.5, 3.0});
  const Axis b({0.0, 1.0, 4.0});

  const MultiAxis ma(a, b);

  // test general properties
  BOOST_CHECK_EQUAL(ma.getNAxes(), 2u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(0), 2u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(1), 2u);
  BOOST_CHECK_EQUAL(ma.getNTotalBins(), 4u);

  // test grid points
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0}), 5u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 1}), 6u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 4}), 7u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.5, 0}), 9u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.5, 1}), 10u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.5, 4}), 11u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3, 0}), 13u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3, 1}), 14u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3, 4}), 15u);

  // test some arbitrary points
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.3, 1.2}), 6u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3.3, 2.2}), 14u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1.8, 0.9}), 9u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({3.1, 0.7}), 13u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2.3, 1.4}), 10u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({2, -3}), 8u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 8}), 11u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-3, 1}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({11, 3}), 14u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-3, -2}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({7, -2}), 12u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-1, 12}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({11, 12}), 15u);

  // flat bin index -> multi bin indices
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(0) == (MultiIndex{0, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(1) == (MultiIndex{0, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(2) == (MultiIndex{0, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(3) == (MultiIndex{0, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(4) == (MultiIndex{1, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(5) == (MultiIndex{1, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(6) == (MultiIndex{1, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(7) == (MultiIndex{1, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(8) == (MultiIndex{2, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(9) == (MultiIndex{2, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(10) == (MultiIndex{2, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(11) == (MultiIndex{2, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(12) == (MultiIndex{3, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(13) == (MultiIndex{3, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(14) == (MultiIndex{3, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(15) == (MultiIndex{3, 3}));

  // multi bin indices -> flat bin index
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 0}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 1}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 2}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 3}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 0}), 4u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 1}), 5u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 2}), 6u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 3}), 7u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 0}), 8u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 1}), 9u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 2}), 10u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 3}), 11u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 0}), 12u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 1}), 13u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 2}), 14u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 3}), 15u);

  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(
                  ma.getGlobalBinFromPoint({3.2, 1.8})) == (MultiIndex{3, 2}));

  // inside checks
  BOOST_CHECK(!ma.isInside({-2., -1}));
  BOOST_CHECK(!ma.isInside({-2., 1.}));
  BOOST_CHECK(!ma.isInside({-2., 5.}));
  BOOST_CHECK(!ma.isInside({1., -1.}));
  BOOST_CHECK(!ma.isInside({6., -1.}));
  BOOST_CHECK(ma.isInside({0.5, 1.3}));
  BOOST_CHECK(!ma.isInside({3., -1.}));
  BOOST_CHECK(!ma.isInside({3., 0.3}));
  BOOST_CHECK(!ma.isInside({3., 4.}));
  BOOST_CHECK(!ma.isInside({-1., 4.}));
  BOOST_CHECK(!ma.isInside({2., 4.}));
  BOOST_CHECK(!ma.isInside({5., 4.}));

  // test some bin centers
  CHECK_CLOSE_ABS(ma.getBinCenter({{1, 1}}), (Point{0.25, 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({{2, 1}}), (Point{1.75, 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({{1, 2}}), (Point{0.25, 2.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({{2, 2}}), (Point{1.75, 2.5}), 1e-6);

  // test some lower-left bin edges
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({{1, 1}}), (Point{0., 0.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({{2, 1}}), (Point{0.5, 0.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({{1, 2}}), (Point{0., 1.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({{2, 2}}), (Point{0.5, 1.}), 1e-6);

  // test some upper right-bin edges
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({{1, 1}}), (Point{0.5, 1.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({{2, 1}}), (Point{3., 1.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({{1, 2}}), (Point{0.5, 4.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({{2, 2}}), (Point{3., 4.}), 1e-6);
}

BOOST_AUTO_TEST_CASE(test_3d_variable) {
  using Point = std::array<double, 3>;
  using MultiIndex = std::array<std::size_t, 3>;

  const Axis a({0.0, 1.0});
  const Axis b({0.0, 0.5, 3.0});
  const Axis c({0.0, 0.5, 3.0, 3.3});

  const MultiAxis ma(a, b, c);

  // test general properties
  BOOST_CHECK_EQUAL(ma.getNAxes(), 3u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(0), 1u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(1), 2u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(2), 3u);
  BOOST_CHECK_EQUAL(ma.getNTotalBins(), 6u);

  // test grid points
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0, 0}), 26u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0, 0}), 46u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0.5, 0}), 31u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0.5, 0}), 51u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 3, 0}), 36u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 3, 0}), 56u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0, 0.5}), 27u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0, 0.5}), 47u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0.5, 0.5}), 32u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0.5, 0.5}), 52u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 3, 0.5}), 37u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 3, 0.5}), 57u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0, 3}), 28u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0, 3}), 48u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0.5, 3}), 33u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0.5, 3}), 53u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 3, 3}), 38u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 3, 3}), 58u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0, 3.3}), 29u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0, 3.3}), 49u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0.5, 3.3}), 34u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0.5, 3.3}), 54u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 3, 3.3}), 39u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 3, 3.3}), 59u);

  // flat bin index -> multi bin indices
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(0) == (MultiIndex{0, 0, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(1) == (MultiIndex{0, 0, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(2) == (MultiIndex{0, 0, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(3) == (MultiIndex{0, 0, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(4) == (MultiIndex{0, 0, 4}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(5) == (MultiIndex{0, 1, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(21) == (MultiIndex{1, 0, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(22) == (MultiIndex{1, 0, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(23) == (MultiIndex{1, 0, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(24) == (MultiIndex{1, 0, 4}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(25) == (MultiIndex{1, 1, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(26) == (MultiIndex{1, 1, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(57) == (MultiIndex{2, 3, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(58) == (MultiIndex{2, 3, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(59) == (MultiIndex{2, 3, 4}));

  // multi bin indices -> flat bin index
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 0, 0}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 0, 0}), 20u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 0, 0}), 40u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 1, 0}), 5u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 1, 0}), 25u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 1, 0}), 45u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 3, 1}), 16u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 3, 1}), 36u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 3, 1}), 56u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 0, 2}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 0, 2}), 22u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 0, 2}), 42u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 3, 4}), 19u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 3, 4}), 39u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 3, 4}), 59u);

  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(ma.getGlobalBinFromPoint(
                  {1.8, 0.7, 3.2})) == (MultiIndex{2, 2, 3}));

  // inside checks
  BOOST_CHECK(!ma.isInside({-2., -1, -2}));
  BOOST_CHECK(!ma.isInside({-2., 1., 0.}));
  BOOST_CHECK(!ma.isInside({-2., 5., -1}));
  BOOST_CHECK(!ma.isInside({1., -1., 1.}));
  BOOST_CHECK(!ma.isInside({6., -1., 4.}));
  BOOST_CHECK(ma.isInside({0.5, 1.3, 1.7}));
  BOOST_CHECK(!ma.isInside({1., -1., -0.4}));
  BOOST_CHECK(!ma.isInside({1., 0.3, 3.4}));
  BOOST_CHECK(!ma.isInside({1., 3., 0.8}));
  BOOST_CHECK(!ma.isInside({-1., 3., 5.}));
  BOOST_CHECK(!ma.isInside({2., 3., -1.}));
  BOOST_CHECK(!ma.isInside({5., 3., 0.5}));

  // test some bin centers
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 1, 1}), (Point{0.5, 0.25, 0.25}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 1, 2}), (Point{0.5, 0.25, 1.75}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 1, 3}), (Point{0.5, 0.25, 3.15}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 2, 1}), (Point{0.5, 1.75, 0.25}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 2, 2}), (Point{0.5, 1.75, 1.75}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 2, 3}), (Point{0.5, 1.75, 3.15}), 1e-6);

  // test some lower-left bin edges
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 1, 1}), (Point{0., 0., 0.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 1, 2}), (Point{0., 0., 0.5}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 1, 3}), (Point{0., 0., 3.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 2, 1}), (Point{0., 0.5, 0.}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 2, 2}), (Point{0., 0.5, 0.5}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 2, 3}), (Point{0., 0.5, 3.}),
                  1e-6);

  // test some upper right-bin edges
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1, 1, 1}), (Point{1., 0.5, 0.5}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1, 1, 2}), (Point{1., 0.5, 3.}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1, 1, 3}), (Point{1., 0.5, 3.3}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1, 2, 1}), (Point{1., 3., 0.5}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1, 2, 2}), (Point{1., 3., 3.}),
                  1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1, 2, 3}), (Point{1., 3., 3.3}),
                  1e-6);
}

BOOST_AUTO_TEST_CASE(test_2d_mixed) {
  using Point = std::array<double, 2>;
  using MultiIndex = std::array<std::size_t, 2>;

  const Axis a(0.0, 1.0, 4u);
  const Axis b({0.0, 0.5, 3.0});

  const MultiAxis ma(a, b);

  // test general properties
  BOOST_CHECK_EQUAL(ma.getNAxes(), 2u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(0), 4u);
  BOOST_CHECK_EQUAL(ma.getNBins().at(1), 2u);
  BOOST_CHECK_EQUAL(ma.getNTotalBins(), 8u);

  // test grid points
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0}), 5u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.25, 0}), 9u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.5, 0}), 13u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.75, 0}), 17u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0}), 21u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 0.5}), 6u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.25, 0.5}), 10u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.5, 0.5}), 14u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.75, 0.5}), 18u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 0.5}), 22u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0, 3}), 7u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.25, 3}), 11u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.5, 3}), 15u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.75, 3}), 19u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1, 3}), 23u);

  // test some arbitrary points
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({1.2, 0.3}), 21u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.2, 1.3}), 6u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.9, 1.8}), 18u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.7, 2.1}), 14u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.4, 0.3}), 9u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-3, 2}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({8, 1}), 22u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.1, -3}), 4u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({0.8, 11}), 19u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-2, -3}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({-2, 7}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({12, -1}), 20u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromPoint({12, 11}), 23u);

  // flat bin index -> multi bin indices
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(0) == (MultiIndex{0, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(1) == (MultiIndex{0, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(2) == (MultiIndex{0, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(3) == (MultiIndex{0, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(4) == (MultiIndex{1, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(5) == (MultiIndex{1, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(6) == (MultiIndex{1, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(7) == (MultiIndex{1, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(8) == (MultiIndex{2, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(9) == (MultiIndex{2, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(10) == (MultiIndex{2, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(11) == (MultiIndex{2, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(12) == (MultiIndex{3, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(13) == (MultiIndex{3, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(14) == (MultiIndex{3, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(15) == (MultiIndex{3, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(16) == (MultiIndex{4, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(17) == (MultiIndex{4, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(18) == (MultiIndex{4, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(19) == (MultiIndex{4, 3}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(20) == (MultiIndex{5, 0}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(21) == (MultiIndex{5, 1}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(22) == (MultiIndex{5, 2}));
  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(23) == (MultiIndex{5, 3}));

  // multi bin indices -> flat bin index
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 0}), 0u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 1}), 1u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 2}), 2u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({0, 3}), 3u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 0}), 4u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 1}), 5u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 2}), 6u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({1, 3}), 7u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 0}), 8u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 1}), 9u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 2}), 10u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({2, 3}), 11u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 0}), 12u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 1}), 13u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 2}), 14u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({3, 3}), 15u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({4, 0}), 16u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({4, 1}), 17u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({4, 2}), 18u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({4, 3}), 19u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({5, 0}), 20u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({5, 1}), 21u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({5, 2}), 22u);
  BOOST_CHECK_EQUAL(ma.getGlobalBinFromLocalBins({5, 3}), 23u);

  BOOST_CHECK(ma.getLocalBinsFromGlobalBin(ma.getGlobalBinFromPoint(
                  Point({{1.1, 1.7}}))) == MultiIndex({{5, 2}}));

  // inside checks
  BOOST_CHECK(!ma.isInside({-2., -1}));
  BOOST_CHECK(!ma.isInside({-2., 1.}));
  BOOST_CHECK(!ma.isInside({-2., 5.}));
  BOOST_CHECK(!ma.isInside({0.1, -1.}));
  BOOST_CHECK(!ma.isInside({6., -1.}));
  BOOST_CHECK(ma.isInside({0.5, 1.3}));
  BOOST_CHECK(!ma.isInside({1., -1.}));
  BOOST_CHECK(!ma.isInside({1., 0.3}));
  BOOST_CHECK(!ma.isInside({1., 3.}));
  BOOST_CHECK(!ma.isInside({-1., 3.}));
  BOOST_CHECK(!ma.isInside({0.2, 3.}));
  BOOST_CHECK(!ma.isInside({5., 3.}));

  // test some bin centers
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 1}), (Point{0.125, 0.25}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({1, 2}), (Point{0.125, 1.75}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({2, 1}), (Point{0.375, 0.25}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({2, 2}), (Point{0.375, 1.75}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({3, 1}), (Point{0.625, 0.25}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({3, 2}), (Point{0.625, 1.75}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({4, 1}), (Point{0.875, 0.25}), 1e-6);
  CHECK_CLOSE_ABS(ma.getBinCenter({4, 2}), (Point{0.875, 1.75}), 1e-6);

  // test some lower-left bin edges
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 1}), (Point{0., 0.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({1, 2}), (Point{0., 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({2, 1}), (Point{0.25, 0.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({2, 2}), (Point{0.25, 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({3, 1}), (Point{0.5, 0.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({3, 2}), (Point{0.5, 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({4, 1}), (Point{0.75, 0.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getLowerLeftBinEdge({4, 2}), (Point{0.75, 0.5}), 1e-6);

  // test some upper-right bin edges
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1, 1}), (Point{0.25, 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({1, 2}), (Point{0.25, 3.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({2, 1}), (Point{0.5, 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({2, 2}), (Point{0.5, 3.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({3, 1}), (Point{0.75, 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({3, 2}), (Point{0.75, 3.}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({4, 1}), (Point{1., 0.5}), 1e-6);
  CHECK_CLOSE_ABS(ma.getUpperRightBinEdge({4, 2}), (Point{1., 3.}), 1e-6);
}

BOOST_AUTO_TEST_CASE(neighborhood) {
  using bins_t = std::vector<std::size_t>;

  const Axis a(0.0, 1.0, 10u);
  const Axis b(0.0, 1.0, 10u);
  const Axis c(0.0, 1.0, 10u);

  const MultiAxis ma1(a);
  const MultiAxis ma2(a, b);
  const MultiAxis ma3(a, b, c);

  // 1D case
  BOOST_CHECK(ma1.getNeighborHoodIndices({0}, 1).collectVector() ==
              (bins_t{0, 1}));
  BOOST_CHECK(ma1.getNeighborHoodIndices({0}, 2).collectVector() ==
              (bins_t{0, 1, 2}));
  BOOST_CHECK(ma1.getNeighborHoodIndices({1}, 1).collectVector() ==
              (bins_t{0, 1, 2}));
  BOOST_CHECK(ma1.getNeighborHoodIndices({1}, 3).collectVector() ==
              (bins_t{0, 1, 2, 3, 4}));
  BOOST_CHECK(ma1.getNeighborHoodIndices({4}, 2).collectVector() ==
              (bins_t{2, 3, 4, 5, 6}));
  BOOST_CHECK(ma1.getNeighborHoodIndices({9}, 2).collectVector() ==
              (bins_t{7, 8, 9, 10, 11}));
  BOOST_CHECK(ma1.getNeighborHoodIndices({10}, 2).collectVector() ==
              (bins_t{8, 9, 10, 11}));
  BOOST_CHECK(ma1.getNeighborHoodIndices({11}, 2).collectVector() ==
              (bins_t{9, 10, 11}));

  // 2D case
  BOOST_CHECK(ma2.getNeighborHoodIndices({0, 0}, 1).collectVector() ==
              (bins_t{0, 1, 12, 13}));
  BOOST_CHECK(ma2.getNeighborHoodIndices({0, 1}, 1).collectVector() ==
              (bins_t{0, 1, 2, 12, 13, 14}));
  BOOST_CHECK(ma2.getNeighborHoodIndices({1, 0}, 1).collectVector() ==
              (bins_t{0, 1, 12, 13, 24, 25}));
  BOOST_CHECK(ma2.getNeighborHoodIndices({1, 1}, 1).collectVector() ==
              (bins_t{0, 1, 2, 12, 13, 14, 24, 25, 26}));
  BOOST_CHECK(ma2.getNeighborHoodIndices({5, 5}, 1).collectVector() ==
              (bins_t{52, 53, 54, 64, 65, 66, 76, 77, 78}));
  BOOST_CHECK(ma2.getNeighborHoodIndices({9, 10}, 2).collectVector() ==
              (bins_t{92,  93,  94,  95,  104, 105, 106, 107, 116, 117,
                      118, 119, 128, 129, 130, 131, 140, 141, 142, 143}));

  // 3D case
  BOOST_CHECK(ma3.getNeighborHoodIndices({0, 0, 0}, 1).collectVector() ==
              (bins_t{0, 1, 12, 13, 144, 145, 156, 157}));
  BOOST_CHECK(ma3.getNeighborHoodIndices({0, 0, 1}, 1).collectVector() ==
              (bins_t{0, 1, 2, 12, 13, 14, 144, 145, 146, 156, 157, 158}));
  BOOST_CHECK(ma3.getNeighborHoodIndices({0, 1, 0}, 1).collectVector() ==
              (bins_t{0, 1, 12, 13, 24, 25, 144, 145, 156, 157, 168, 169}));
  BOOST_CHECK(ma3.getNeighborHoodIndices({1, 0, 0}, 1).collectVector() ==
              (bins_t{0, 1, 12, 13, 144, 145, 156, 157, 288, 289, 300, 301}));
  BOOST_CHECK(ma3.getNeighborHoodIndices({0, 1, 1}, 1).collectVector() ==
              (bins_t{0, 1, 2, 12, 13, 14, 24, 25, 26, 144, 145, 146, 156, 157,
                      158, 168, 169, 170}));
  BOOST_CHECK(ma3.getNeighborHoodIndices({1, 1, 1}, 1).collectVector() ==
              (bins_t{0,   1,   2,   12,  13,  14,  24,  25,  26,
                      144, 145, 146, 156, 157, 158, 168, 169, 170,
                      288, 289, 290, 300, 301, 302, 312, 313, 314}));
  BOOST_CHECK(ma3.getNeighborHoodIndices({11, 10, 9}, 1).collectVector() ==
              (bins_t{1556, 1557, 1558, 1568, 1569, 1570, 1580, 1581, 1582,
                      1700, 1701, 1702, 1712, 1713, 1714, 1724, 1725, 1726}));

  // Neighbors array
  std::array<std::pair<int, int>, 1> a1;
  a1.at(0) = std::make_pair<int, int>(-1, 1);
  BOOST_CHECK(ma1.getNeighborHoodIndices({0}, a1).collectVector() ==
              (bins_t{0, 1}));
  BOOST_CHECK(ma1.getNeighborHoodIndices({2}, a1).collectVector() ==
              (bins_t{1, 2, 3}));

  a1.at(0) = std::make_pair<int, int>(2, 3);
  BOOST_CHECK(ma1.getNeighborHoodIndices({2}, a1).collectVector() ==
              (bins_t{4, 5}));

  a1.at(0) = std::make_pair<int, int>(-2, -1);
  BOOST_CHECK(ma1.getNeighborHoodIndices({2}, a1).collectVector() ==
              (bins_t{0, 1}));

  a1.at(0) = std::make_pair<int, int>(-3, -1);
  BOOST_CHECK(ma1.getNeighborHoodIndices({2}, a1).collectVector() ==
              (bins_t{0, 1}));

  const Axis d(AxisClosed, 0.0, 1.0, 10u);

  const MultiAxis ma1Cl(d);

  BOOST_CHECK(ma1Cl.getNeighborHoodIndices({0}, 1)
                  .collectVector()
                  .empty());  // underflow, makes no sense
  BOOST_CHECK(ma1Cl.getNeighborHoodIndices({11}, 1)
                  .collectVector()
                  .empty());  // overflow, makes no sense
  BOOST_CHECK(ma1Cl.getNeighborHoodIndices({1}, 1).collectVector() ==
              (bins_t{9, 0, 1}));
  BOOST_CHECK(ma1Cl.getNeighborHoodIndices({5}, 1).collectVector() ==
              (bins_t{3, 4, 5}));

  const Axis f(AxisClosed, 0.0, 1.0, 5u);
  const Axis e(AxisClosed, 0.0, 1.0, 5u);

  const MultiAxis ma2Cl(e, f);

  BOOST_CHECK(ma2Cl.getNeighborHoodIndices({3, 3}, 1).collectVector() ==
              (bins_t{6, 7, 8, 11, 12, 13, 16, 17, 18}));
  BOOST_CHECK(ma2Cl.getNeighborHoodIndices({1, 1}, 1).collectVector() ==
              (bins_t{24, 20, 21, 4, 0, 1, 9, 5, 6}));
  BOOST_CHECK(ma2Cl.getNeighborHoodIndices({1, 5}, 1).collectVector() ==
              (bins_t{23, 24, 20, 3, 4, 0, 8, 9, 5}));
  BOOST_CHECK(ma2Cl.getNeighborHoodIndices({5, 1}, 1).collectVector() ==
              (bins_t{19, 15, 16, 24, 20, 21, 4, 0, 1}));
  BOOST_CHECK(ma2Cl.getNeighborHoodIndices({5, 5}, 1).collectVector() ==
              (bins_t{18, 19, 15, 23, 24, 20, 3, 4, 0}));

  BOOST_CHECK(ma2Cl.getNeighborHoodIndices({3, 3}, 2).collectVector() ==
              (bins_t{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                      13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24}));
  BOOST_CHECK(ma2Cl.getNeighborHoodIndices({1, 1}, 2).collectVector() ==
              (bins_t{18, 19, 15, 16, 17, 23, 24, 20, 21, 22, 3,  4, 0,
                      1,  2,  8,  9,  5,  6,  7,  13, 14, 10, 11, 12}));
  BOOST_CHECK(ma2Cl.getNeighborHoodIndices({1, 5}, 2).collectVector() ==
              (bins_t{17, 18, 19, 15, 16, 22, 23, 24, 20, 21, 2,  3, 4,
                      0,  1,  7,  8,  9,  5,  6,  12, 13, 14, 10, 11}));
  BOOST_CHECK(ma2Cl.getNeighborHoodIndices({5, 1}, 2).collectVector() ==
              (bins_t{13, 14, 10, 11, 12, 18, 19, 15, 16, 17, 23, 24, 20,
                      21, 22, 3,  4,  0,  1,  2,  8,  9,  5,  6,  7}));
  BOOST_CHECK(ma2Cl.getNeighborHoodIndices({5, 5}, 2).collectVector() ==
              (bins_t{12, 13, 14, 10, 11, 17, 18, 19, 15, 16, 22, 23, 24,
                      20, 21, 2,  3,  4,  0,  1,  7,  8,  9,  5,  6}));

  std::array<std::pair<int, int>, 2> a2;
  a2.at(0) =
      std::make_pair<int, int>(-2, -1);  // only 2 bins left of requested bin
                                         // (not including the requested bin)
  a2.at(1) = std::make_pair<int, int>(
      -1, 2);  // one bin left of requested bin, the requested bin itself and 2
               // bins right of requested bin
  std::set<std::size_t> returnedBins;

  auto returnedBinsVec =
      ma2Cl.getNeighborHoodIndices({3, 2}, a2).collectVector();
  returnedBins.insert(returnedBinsVec.begin(), returnedBinsVec.end());
  std::set<std::size_t> expectedBins{{0, 1, 2, 3, 5, 6, 7, 8}};
  BOOST_CHECK(returnedBins == expectedBins);

  returnedBinsVec = ma2Cl.getNeighborHoodIndices({1, 5}, a2).collectVector();
  returnedBins.clear();
  returnedBins.insert(returnedBinsVec.begin(), returnedBinsVec.end());
  expectedBins = {{15, 16, 18, 19, 20, 21, 23, 24}};
  BOOST_CHECK(returnedBins == expectedBins);

  a2.at(0) = {-6, 7};
  a2.at(1) = {0, 0};
  returnedBinsVec = ma2Cl.getNeighborHoodIndices({1, 5}, a2).collectVector();
  returnedBins.clear();
  returnedBins.insert(returnedBinsVec.begin(), returnedBinsVec.end());
  expectedBins = {{4, 9, 14, 19, 24}};
  BOOST_CHECK(returnedBins == expectedBins);

  // @TODO 3D test would be nice, but should essentially not be a problem if
  // 2D works.

  // clang-format off
  /*
   *       1   2    3    4    5
   *   |------------------------|
   * 1 |  0 |  1 |  2 |  3 |  4 |
   *   |----|----|----|----|----|
   * 2 |  5 |  6 |  7 |  8 |  9 |
   *   |----|----|----|----|----|
   * 3 | 10 | 11 | 12 | 13 | 14 |
   *   |----|----|----|----|----|
   * 4 | 15 | 16 | 17 | 18 | 19 |
   *   |----|----|----|----|----|
   * 5 | 20 | 21 | 22 | 23 | 24 |
   *   |------------------------|
   */
  // clang-format on
}

BOOST_AUTO_TEST_CASE(closestPoints) {
  using Point1 = std::array<double, 1>;
  using Point2 = std::array<double, 2>;
  using Point3 = std::array<double, 3>;
  using bins_t = std::vector<std::size_t>;

  const Axis a(0.0, 1.0, 10u);
  const Axis b(0.0, 1.0, 5u);
  const Axis c(0.0, 1.0, 3u);

  const MultiAxis ma1(a);
  const MultiAxis ma2(a, b);
  const MultiAxis ma3(a, b, c);

  // 1D case
  CHECK_EQUAL_COLLECTIONS(
      ma1.getClosestPointsIndices(Point1{0.52}).collectVector(),
      (bins_t{6, 7}));
  CHECK_EQUAL_COLLECTIONS(
      ma1.getClosestPointsIndices(Point1{0.98}).collectVector(),
      (bins_t{10, 11}));

  // 2D case
  CHECK_EQUAL_COLLECTIONS(
      ma2.getClosestPointsIndices(Point2{0.52, 0.08}).collectVector(),
      (bins_t{43, 44, 50, 51}));
  CHECK_EQUAL_COLLECTIONS(
      ma2.getNeighborHoodIndices(ma2.getLocalBinsFromPoint(Point2{0.05, 0.08}),
                                 {0, 1})
          .collectVector(),
      (bins_t{8, 9, 15, 16}));

  // 3D case
  CHECK_EQUAL_COLLECTIONS(
      ma3.getClosestPointsIndices(Point3{0.23, 0.13, 0.61}).collectVector(),
      (bins_t{112, 113, 117, 118, 147, 148, 152, 153}));
  CHECK_EQUAL_COLLECTIONS(
      ma3.getClosestPointsIndices(Point3{0.52, 0.35, 0.71}).collectVector(),
      (bins_t{223, 224, 228, 229, 258, 259, 263, 264}));

  using EAxisClosed = Axis<AxisType::Equidistant, AxisBoundaryType::Closed>;

  using MultiAxis1Cl_t = MultiAxis<EAxisClosed>;
  using MultiAxis2Cl_t = MultiAxis<EAxisClosed, EAxisClosed>;
  // using MultiAxis3Cl_t = MultiAxis<EAxisClosed, EAxisClosed, EAxisClosed>;

  EAxisClosed aCl(0.0, 1.0, 10u);
  EAxisClosed bCl(0.0, 1.0, 5u);
  // EAxisClosed   cCl(0.0, 1.0, 3u);

  MultiAxis1Cl_t ma1Cl(std::make_tuple(aCl));
  MultiAxis2Cl_t ma2Cl(std::make_tuple(aCl, bCl));

  // 1D case
  CHECK_EQUAL_COLLECTIONS(
      ma1Cl.getClosestPointsIndices(Point1{0.52}).collectVector(),
      (bins_t{5, 6}));
  CHECK_EQUAL_COLLECTIONS(
      ma1Cl.getClosestPointsIndices(Point1{0.98}).collectVector(),
      (bins_t{9, 0}));

  // 2D case
  CHECK_EQUAL_COLLECTIONS(
      ma2Cl.getClosestPointsIndices(Point2{0.52, 0.08}).collectVector(),
      (bins_t{25, 26, 30, 31}));
  CHECK_EQUAL_COLLECTIONS(
      ma2Cl.getClosestPointsIndices(Point2{0.52, 0.68}).collectVector(),
      (bins_t{28, 29, 33, 34}));
  CHECK_EQUAL_COLLECTIONS(
      ma2Cl.getClosestPointsIndices(Point2{0.52, 0.88}).collectVector(),
      (bins_t{29, 25, 34, 30}));
  CHECK_EQUAL_COLLECTIONS(
      ma2Cl.getClosestPointsIndices(Point2{0.05, 0.08}).collectVector(),
      (bins_t{0, 1, 5, 6}));
  CHECK_EQUAL_COLLECTIONS(
      ma2Cl.getClosestPointsIndices(Point2{0.9, 0.95}).collectVector(),
      (bins_t{49, 45, 4, 0}));

  // @TODO: 3D checks would also be nice

  const Axis aOp(AxisBound, 0.0, 1.0, 10u);
  const Axis bOp(AxisBound, 0.0, 1.0, 5u);

  const MultiAxis ma1Op(aOp);
  const MultiAxis ma2Op(aOp, bOp);

  // 1D case
  CHECK_EQUAL_COLLECTIONS(
      ma1Op.getClosestPointsIndices(Point1{0.52}).collectVector(),
      (bins_t{5, 6}));
  CHECK_EQUAL_COLLECTIONS(
      ma1Op.getClosestPointsIndices(Point1{0.98}).collectVector(), (bins_t{9}));
  CHECK_EQUAL_COLLECTIONS(
      ma1Op.getClosestPointsIndices(Point1{0.88}).collectVector(),
      (bins_t{8, 9}));

  // 2D case
  CHECK_EQUAL_COLLECTIONS(
      ma2Op.getClosestPointsIndices(Point2{0.52, 0.08}).collectVector(),
      (bins_t{25, 26, 30, 31}));
  CHECK_EQUAL_COLLECTIONS(
      ma2Op.getClosestPointsIndices(Point2{0.52, 0.68}).collectVector(),
      (bins_t{28, 29, 33, 34}));
  CHECK_EQUAL_COLLECTIONS(
      ma2Op.getClosestPointsIndices(Point2{0.52, 0.88}).collectVector(),
      (bins_t{29, 34}));
  CHECK_EQUAL_COLLECTIONS(
      ma2Op.getClosestPointsIndices(Point2{0.05, 0.1}).collectVector(),
      (bins_t{0, 1, 5, 6}));
  CHECK_EQUAL_COLLECTIONS(
      ma2Op.getClosestPointsIndices(Point2{0.95, 0.95}).collectVector(),
      (bins_t{49}));

  // @TODO: 3D checks would also be nice

  // clang-format off
  /*
   *       1    2    3    4    5
   *    |------------------------|
   *  1 |  0 |  1 |  2 |  3 |  4 |
   *    |----|----|----|----|----|
   *  2 |  5 |  6 |  7 |  8 |  9 |
   *    |----|----|----|----|----|
   *  3 | 10 | 11 | 12 | 13 | 14 |
   *    |----|----|----|----|----|
   *  4 | 15 | 16 | 17 | 18 | 19 |
   *    |----|----|----|----|----|
   *  5 | 20 | 21 | 22 | 23 | 24 |
   *    |------------------------|
   *  6 | 25 | 26 | 27 | 28 | 29 |
   *    |------------------------|
   *  7 | 30 | 31 | 32 | 33 | 34 |
   *    |------------------------|
   *  8 | 35 | 36 | 37 | 38 | 39 |
   *    |------------------------|
   *  9 | 40 | 41 | 42 | 43 | 44 |
   *    |------------------------|
   * 10 | 45 | 46 | 47 | 48 | 49 |
   *    |------------------------|
   */
  // clang-format on
}

BOOST_AUTO_TEST_CASE(Output) {
  const Axis a{AxisOpen, 0.0, 1.0, 10u};
  const Axis b{AxisBound, {1, 2, 3}};

  const MultiAxis ma(a, b);

  std::stringstream ss;
  ss << ma;
  BOOST_CHECK_EQUAL(ss.str(),
                    "Axis<Equidistant, Open>(0, 1, 10, Undefined), "
                    "Axis<Variable, Bound>({1, 2, 3}, Undefined)");

  const IMultiAxis& ima = ma;

  ss.str("");

  ss << ima;

  BOOST_CHECK_EQUAL(ss.str(),
                    "Axis<Equidistant, Open>(0, 1, 10, Undefined), "
                    "Axis<Variable, Bound>({1, 2, 3}, Undefined)");
}

BOOST_AUTO_TEST_CASE(Equality) {
  const Axis a{AxisOpen, 0.0, 1.0, 10u};
  const Axis b{AxisBound, {1, 2, 3}};
  const Axis c{AxisClosed, {1, 2, 5}};

  const MultiAxis ma_ab(a, b);
  const MultiAxis ma_ac(a, c);

  BOOST_CHECK_EQUAL(ma_ab, ma_ab);
  BOOST_CHECK_EQUAL(ma_ac, ma_ac);
  BOOST_CHECK_NE(ma_ab, ma_ac);

  const IMultiAxis& ima_ab = ma_ab;
  const IMultiAxis& ima_ac = ma_ac;

  BOOST_CHECK_EQUAL(ima_ab, ima_ab);
  BOOST_CHECK_EQUAL(ima_ac, ima_ac);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
