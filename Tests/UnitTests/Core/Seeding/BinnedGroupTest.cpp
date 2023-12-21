// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/BinnedGroup.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/GridIterator.hpp"

#include <array>
#include <memory>
#include <vector>

namespace Acts::Test {

BOOST_AUTO_TEST_CASE(binned_group_constructor) {
  Acts::detail::EquidistantAxis xAxis(0, 100, 10);
  Acts::detail::EquidistantAxis yAxis(0, 100, 10);
  Acts::detail::EquidistantAxis zAxis(0, 100, 10);

  using grid_1d_t =
      Acts::Grid<std::vector<std::size_t>, Acts::detail::EquidistantAxis>;
  using grid_2d_t =
      Acts::Grid<std::vector<std::size_t>, Acts::detail::EquidistantAxis,
                 Acts::detail::EquidistantAxis>;
  using grid_3d_t =
      Acts::Grid<std::vector<std::size_t>, Acts::detail::EquidistantAxis,
                 Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis>;

  std::unique_ptr<grid_1d_t> grid_1d =
      std::make_unique<grid_1d_t>(std::make_tuple(xAxis));
  std::unique_ptr<grid_2d_t> grid_2d =
      std::make_unique<grid_2d_t>(std::make_tuple(xAxis, yAxis));
  std::unique_ptr<grid_3d_t> grid_3d =
      std::make_unique<grid_3d_t>(std::make_tuple(xAxis, yAxis, zAxis));

  std::shared_ptr<Acts::GridBinFinder<1ul>> binFinder_1d =
      std::make_shared<Acts::GridBinFinder<1ul>>(1);
  std::shared_ptr<Acts::GridBinFinder<2ul>> binFinder_2d =
      std::make_shared<Acts::GridBinFinder<2ul>>(1, 1);
  std::shared_ptr<Acts::GridBinFinder<3ul>> binFinder_3d =
      std::make_shared<Acts::GridBinFinder<3ul>>(1, 2, 1);

  std::array<std::vector<std::size_t>, 1ul> navigation_1d;
  navigation_1d[0ul].resize(10);
  std::array<std::vector<std::size_t>, 2ul> navigation_2d;
  navigation_2d[1ul].resize(10);

  // Costructors
  // We provide a proper navigation
  Acts::BinnedGroup<grid_1d_t> group_1d(std::move(grid_1d), binFinder_1d,
                                        binFinder_1d, std::move(navigation_1d));
  // We provide a partial navigation, the constructor will complete it
  Acts::BinnedGroup<grid_2d_t> group_2d(std::move(grid_2d), binFinder_2d,
                                        binFinder_2d, std::move(navigation_2d));
  // We do not provide navigation, the constructor will define it
  Acts::BinnedGroup<grid_3d_t> group_3d(std::move(grid_3d), binFinder_3d,
                                        binFinder_3d);

  BOOST_CHECK_NO_THROW(group_1d.grid());
  BOOST_CHECK_NO_THROW(group_2d.grid());
  BOOST_CHECK_NO_THROW(group_3d.grid());

  // Move Constructor/Assignment
  const Acts::BinnedGroup<grid_1d_t> group_1d_moved(std::move(group_1d));
  const Acts::BinnedGroup<grid_2d_t> group_2d_moved(std::move(group_2d));
  const Acts::BinnedGroup<grid_3d_t> group_3d_moved = std::move(group_3d);

  BOOST_CHECK_THROW(group_1d.grid(), std::runtime_error);
  BOOST_CHECK_THROW(group_2d.grid(), std::runtime_error);
  BOOST_CHECK_THROW(group_3d.grid(), std::runtime_error);

  BOOST_CHECK_NO_THROW(group_1d_moved.grid());
  BOOST_CHECK_NO_THROW(group_2d_moved.grid());
  BOOST_CHECK_NO_THROW(group_3d_moved.grid());

  // With some nullptr inputs
  std::unique_ptr<grid_1d_t> grid_1d_nullTest_1 =
      std::make_unique<grid_1d_t>(std::make_tuple(xAxis));
  std::unique_ptr<grid_1d_t> grid_1d_nullTest_2 =
      std::make_unique<grid_1d_t>(std::make_tuple(xAxis));

  Acts::BinnedGroup<grid_1d_t> group_1d_nullBotFinder(
      std::move(grid_1d_nullTest_1), nullptr, binFinder_1d);
  Acts::BinnedGroup<grid_1d_t> group_1d_nullTopFinder(
      std::move(grid_1d_nullTest_2), binFinder_1d, nullptr);

  BOOST_CHECK_NO_THROW(group_1d_nullBotFinder.grid());
  BOOST_CHECK_NO_THROW(group_1d_nullTopFinder.grid());
}

BOOST_AUTO_TEST_CASE(binned_group_iterations_1d_emptyGrid) {
  using grid_t =
      Acts::Grid<std::vector<std::size_t>, Acts::detail::EquidistantAxis>;
  using binfinder_t = Acts::GridBinFinder<1ul>;

  Acts::detail::EquidistantAxis xAxis(0, 100, 10);
  std::unique_ptr<grid_t> grid =
      std::make_unique<grid_t>(std::make_tuple(xAxis));
  std::shared_ptr<binfinder_t> binfinder = std::make_shared<binfinder_t>(0);
  Acts::BinnedGroup<grid_t> group(std::move(grid), binfinder, binfinder);

  Acts::BinnedGroupIterator<grid_t> itrStart = group.begin();
  Acts::BinnedGroupIterator<grid_t> itrStop = group.end();

  std::size_t nIterations = 0ul;
  for (; itrStart != itrStop; ++itrStart, ++nIterations) {
  }

  BOOST_CHECK_EQUAL(nIterations, 0ul);
}

BOOST_AUTO_TEST_CASE(binned_group_iterations_2d_emptyGrid) {
  using grid_t =
      Acts::Grid<std::vector<std::size_t>, Acts::detail::EquidistantAxis,
                 Acts::detail::EquidistantAxis>;
  using binfinder_t = Acts::GridBinFinder<2ul>;

  Acts::detail::EquidistantAxis xAxis(0, 100, 10);
  Acts::detail::EquidistantAxis yAxis(0, 100, 10);
  std::unique_ptr<grid_t> grid = std::make_unique<grid_t>(
      std::make_tuple(std::move(xAxis), std::move(yAxis)));
  std::shared_ptr<binfinder_t> binfinder = std::make_shared<binfinder_t>(0, 0);
  Acts::BinnedGroup<grid_t> group(std::move(grid), binfinder, binfinder);

  Acts::BinnedGroupIterator<grid_t> itrStart = group.begin();
  Acts::BinnedGroupIterator<grid_t> itrStop = group.end();

  std::size_t nIterations = 0ul;
  for (; itrStart != itrStop; ++itrStart, ++nIterations) {
  }

  BOOST_CHECK_EQUAL(nIterations, 0ul);
}

BOOST_AUTO_TEST_CASE(binned_group_iterations_1d_perFilledGrid) {
  using grid_t =
      Acts::Grid<std::vector<std::size_t>, Acts::detail::EquidistantAxis>;
  using binfinder_t = Acts::GridBinFinder<1ul>;

  Acts::detail::EquidistantAxis xAxis(0, 100, 10);
  std::unique_ptr<grid_t> grid =
      std::make_unique<grid_t>(std::make_tuple(xAxis));
  /// Add some entries to the grid
  grid->at(1ul).push_back(4ul);
  grid->at(1ul).push_back(1ul);
  grid->at(8ul).push_back(7ul);
  grid->at(9ul).push_back(2ul);

  std::shared_ptr<binfinder_t> botBinfinder = std::make_shared<binfinder_t>(0);
  std::shared_ptr<binfinder_t> topBinfinder = std::make_shared<binfinder_t>(1);
  Acts::BinnedGroup<grid_t> group(std::move(grid), botBinfinder, topBinfinder);

  std::size_t nIterations = 0ul;
  for (const auto [bottom, middle, top] : group) {
    ++nIterations;
    BOOST_CHECK_EQUAL(bottom.size(), 1ul);
    BOOST_CHECK_EQUAL(top.size(), 3ul);
  }

  BOOST_CHECK_EQUAL(nIterations, 3ul);
}

BOOST_AUTO_TEST_CASE(binned_group_iterations_2d_perFilledGrid) {
  using grid_t =
      Acts::Grid<std::vector<std::size_t>, Acts::detail::EquidistantAxis,
                 Acts::detail::EquidistantAxis>;
  using binfinder_t = Acts::GridBinFinder<2ul>;

  Acts::detail::EquidistantAxis xAxis(0, 100, 10);
  Acts::detail::EquidistantAxis yAxis(0, 100, 10);
  std::unique_ptr<grid_t> grid =
      std::make_unique<grid_t>(std::make_tuple(xAxis, yAxis));
  /// Add some entries to the grid
  grid->atLocalBins({2ul, 4ul}).push_back(4ul);
  grid->atLocalBins({4ul, 4ul}).push_back(4ul);

  std::shared_ptr<binfinder_t> botBinfinder =
      std::make_shared<binfinder_t>(1, 2);
  std::shared_ptr<binfinder_t> topBinfinder =
      std::make_shared<binfinder_t>(1, 1);
  Acts::BinnedGroup<grid_t> group(std::move(grid), botBinfinder, topBinfinder);

  std::size_t nIterations = 0ul;
  for (const auto [bottom, middle, top] : group) {
    ++nIterations;
    BOOST_CHECK_EQUAL(bottom.size(), 15ul);
    BOOST_CHECK_EQUAL(top.size(), 9ul);
  }

  BOOST_CHECK_EQUAL(nIterations, 2ul);
}

BOOST_AUTO_TEST_CASE(binned_group_fill_2d) {
  struct point {};

  using value_t = std::unique_ptr<InternalSpacePoint<point>>;
  using phiAxis_t = Acts::detail::Axis<detail::AxisType::Equidistant,
                                       detail::AxisBoundaryType::Closed>;
  using zAxis_t = detail::Axis<detail::AxisType::Equidistant,
                               detail::AxisBoundaryType::Bound>;
  using grid_t = Acts::Grid<std::vector<value_t>, phiAxis_t, zAxis_t>;
  using binfinder_t = Acts::GridBinFinder<2ul>;

  phiAxis_t phiAxis(-M_PI, M_PI, 40);
  zAxis_t zAxis(0, 100, 10);

  std::unique_ptr<grid_t> grid = std::make_unique<grid_t>(
      std::make_tuple(std::move(phiAxis), std::move(zAxis)));
  std::shared_ptr<binfinder_t> binfinder = std::make_shared<binfinder_t>(1, 1);
  Acts::BinnedGroup<grid_t> group(std::move(grid), binfinder, binfinder);

  auto extractionFunction =
      [](const point&, float, float,
         float) -> std::pair<Acts::Vector3, Acts::Vector2> {
    return std::make_pair(Acts::Vector3({2, 4, 2}), Acts::Vector2({0, 0}));
  };
  Acts::Extent extent;
  Acts::SeedFinderConfig<point> config;
  Acts::SeedFinderOptions options;

  std::vector<const point*> collection;
  for (std::size_t i(0ul); i < 30ul; ++i) {
    collection.push_back(new point());
  }

  BOOST_CHECK_THROW(
      group.fill(config.toInternalUnits(), options.toInternalUnits(),
                 collection.begin(), collection.end(),
                 std::move(extractionFunction), extent),
      std::runtime_error);

  Acts::SeedFilterConfig seedFilterConfig;
  config.seedFilter = std::make_shared<Acts::SeedFilter<point>>(
      seedFilterConfig.toInternalUnits());

  BOOST_CHECK_NO_THROW(group.fill(
      config.toInternalUnits(), options.toInternalUnits(), collection.begin(),
      collection.end(), std::move(extractionFunction), extent));

  BOOST_CHECK_THROW(
      group.fill(config, options.toInternalUnits(), collection.begin(),
                 collection.begin(), std::move(extractionFunction), extent),
      std::runtime_error);

  BOOST_CHECK_THROW(
      group.fill(config.toInternalUnits(), options, collection.begin(),
                 collection.begin(), std::move(extractionFunction), extent),
      std::runtime_error);

  std::size_t nIterations = 0ul;
  for (const auto [bottom, middle, top] : group) {
    ++nIterations;
    const auto& coll = group.grid().at(middle);
    BOOST_CHECK_EQUAL(coll.size(), 30ul);
  }

  BOOST_CHECK_EQUAL(nIterations, 1ul);
}

}  // namespace Acts::Test
