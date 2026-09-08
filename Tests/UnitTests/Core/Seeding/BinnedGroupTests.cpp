// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace Acts::Test {

namespace {

using AxisType = Axis<Acts::AxisType::Equidistant, AxisBoundaryType::Open>;
using GridType = Grid<std::vector<std::size_t>, AxisType, AxisType>;
using BinFinderType = GridBinFinder<2ul>;

// 3 x 4 bins, all filled so the group iterator does not skip any of them
GridType makeFilledGrid() {
  AxisType xAxis(AxisOpen, 0., 3., 3ul);
  AxisType yAxis(AxisOpen, 0., 4., 4ul);
  GridType grid(std::make_tuple(std::move(xAxis), std::move(yAxis)));
  for (std::size_t i(1ul); i <= 3ul; ++i) {
    for (std::size_t j(1ul); j <= 4ul; ++j) {
      grid.atLocalBins({i, j}).push_back(1ul);
    }
  }
  return grid;
}

// Local bin indices of the middle candidates, in visiting order
std::vector<std::array<std::size_t, 2ul>> visitedBins(
    const BinnedGroup<GridType>& group) {
  std::vector<std::array<std::size_t, 2ul>> visited;
  for (const auto& [bottom, middle, top] : group) {
    static_cast<void>(bottom);
    static_cast<void>(top);
    visited.push_back(
        group.grid().multiAxis().getLocalBinsFromGlobalBin(middle));
  }
  return visited;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(BinnedGroupTests)

BOOST_AUTO_TEST_CASE(EmptyNavigationVisitsAllBinsInOrder) {
  BinFinderType bottomFinder(1, 1);
  BinFinderType topFinder(1, 1);
  BinnedGroup<GridType> group(makeFilledGrid(), bottomFinder, topFinder);

  const auto visited = visitedBins(group);
  BOOST_REQUIRE_EQUAL(visited.size(), 12ul);
  BOOST_CHECK_EQUAL(visited.front()[0], 1ul);
  BOOST_CHECK_EQUAL(visited.front()[1], 1ul);
  BOOST_CHECK_EQUAL(visited.back()[0], 3ul);
  BOOST_CHECK_EQUAL(visited.back()[1], 4ul);
}

BOOST_AUTO_TEST_CASE(CustomNavigationVisitsOnlyRequestedBins) {
  BinFinderType bottomFinder(1, 1);
  BinFinderType topFinder(1, 1);
  // Skip the first bin of each axis, reverse the order of the second one
  std::array<std::vector<std::size_t>, 2ul> navigation{
      std::vector<std::size_t>{3ul, 2ul}, std::vector<std::size_t>{4ul, 2ul}};
  BinnedGroup<GridType> group(makeFilledGrid(), bottomFinder, topFinder,
                              navigation);

  const std::vector<std::array<std::size_t, 2ul>> expected{
      {3ul, 4ul}, {3ul, 2ul}, {2ul, 4ul}, {2ul, 2ul}};
  BOOST_CHECK(visitedBins(group) == expected);
}

BOOST_AUTO_TEST_CASE(NavigationRejectsUnderflowBin) {
  BinFinderType bottomFinder(1, 1);
  BinFinderType topFinder(1, 1);
  // 0 is the underflow bin, which is never filled
  std::array<std::vector<std::size_t>, 2ul> navigation{
      std::vector<std::size_t>{0ul, 1ul, 2ul}, std::vector<std::size_t>{}};
  BOOST_CHECK_THROW(BinnedGroup<GridType>(makeFilledGrid(), bottomFinder,
                                          topFinder, navigation),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(NavigationRejectsOverflowBin) {
  BinFinderType bottomFinder(1, 1);
  BinFinderType topFinder(1, 1);
  // The first axis has 3 bins, so 4 is the overflow bin
  std::array<std::vector<std::size_t>, 2ul> navigation{
      std::vector<std::size_t>{1ul, 4ul}, std::vector<std::size_t>{}};
  BOOST_CHECK_THROW(BinnedGroup<GridType>(makeFilledGrid(), bottomFinder,
                                          topFinder, navigation),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(NavigationRejectsDuplicatedBin) {
  BinFinderType bottomFinder(1, 1);
  BinFinderType topFinder(1, 1);
  std::array<std::vector<std::size_t>, 2ul> navigation{
      std::vector<std::size_t>{}, std::vector<std::size_t>{2ul, 3ul, 2ul}};
  BOOST_CHECK_THROW(BinnedGroup<GridType>(makeFilledGrid(), bottomFinder,
                                          topFinder, navigation),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(NavigationIsValidatedWithMaskConstructor) {
  BinFinderType bottomFinder(1, 1);
  BinFinderType topFinder(1, 1);
  GridType grid = makeFilledGrid();
  std::vector<bool> mask(grid.size(true), true);
  std::array<std::vector<std::size_t>, 2ul> navigation{
      std::vector<std::size_t>{0ul}, std::vector<std::size_t>{}};
  BOOST_CHECK_THROW(BinnedGroup<GridType>(std::move(grid), std::move(mask),
                                          bottomFinder, topFinder, navigation),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
