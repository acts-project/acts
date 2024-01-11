// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

using namespace Acts;

namespace bd = boost::unit_test::data;

BOOST_AUTO_TEST_SUITE(GridAccessHelpersTests)

BOOST_AUTO_TEST_CASE(Grid1DAccess) {
  Acts::GridAxisGenerators::EqBound eqBound{{0., 10.}, 10};
  using GridType = Acts::GridAxisGenerators::EqBound::grid_type<std::size_t>;
  using PointType = GridType::point_t;
  auto grid = GridType(eqBound());

  for (std::size_t i = 0; i < 10; ++i) {
    grid.atPosition(PointType{i + 0.5}) = i;
  }

  // Local access
  std::vector<std::size_t> fAccessor = {0u};
  std::vector<std::size_t> sAccessor = {1u};

  Vector2 lPosition{3.5, 6.5};
  auto flAccess =
      GridAccessHelpers::accessLocal<GridType>(lPosition, fAccessor);
  auto slAccess =
      GridAccessHelpers::accessLocal<GridType>(lPosition, sAccessor);

  // This should take out local 1D either first or second
  BOOST_CHECK_EQUAL(grid.atPosition(flAccess), 3u);
  BOOST_CHECK_EQUAL(grid.atPosition(slAccess), 6u);

  // Global access
  Vector3 gPosition{0.5, 3.5, 6.5};
  std::vector<BinningValue> fCast = {Acts::binX};
  std::vector<BinningValue> sCast = {Acts::binY};
  std::vector<BinningValue> tCast = {Acts::binZ};

  auto fgAccess = GridAccessHelpers::castPosition<GridType>(gPosition, fCast);
  auto sgAccess = GridAccessHelpers::castPosition<GridType>(gPosition, sCast);
  auto tgAccess = GridAccessHelpers::castPosition<GridType>(gPosition, tCast);
  BOOST_CHECK_EQUAL(grid.atPosition(fgAccess), 0u);
  BOOST_CHECK_EQUAL(grid.atPosition(sgAccess), 3u);
  BOOST_CHECK_EQUAL(grid.atPosition(tgAccess), 6u);

  // Binned access
  Vector2 flbPosition = GridAccessHelpers::toLocal(grid, 1u, 0u, 0u);
  auto flbAccess =
      GridAccessHelpers::accessLocal<GridType>(flbPosition, fAccessor);
  BOOST_CHECK_EQUAL(grid.atPosition(flbAccess), 1u);
  BOOST_CHECK_THROW(GridAccessHelpers::toLocal(grid, 1u, 0u, 3u),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(Grid2DAccess) {
  Acts::GridAxisGenerators::EqBoundEqBound eqeqBound{
      {0., 10.}, 10, {0., 10.}, 10};
  using GridType =
      Acts::GridAxisGenerators::EqBoundEqBound::grid_type<std::size_t>;
  using PointType = GridType::point_t;
  auto grid = GridType(eqeqBound());
  for (std::size_t j = 0; j < 10u; ++j) {
    for (std::size_t i = 0; i < 10u; ++i) {
      grid.atPosition(PointType{i + 0.5, j + 0.5}) = j * 100 + i;
    }
  }

  // Local access
  std::vector<std::size_t> fAccessor = {0u, 1u};
  Vector2 lPosition{3.5, 6.5};
  auto flAccess =
      GridAccessHelpers::accessLocal<GridType>(lPosition, fAccessor);
  BOOST_CHECK_EQUAL(grid.atPosition(flAccess), 603u);

  // Global access
  Vector3 gPosition{0.5, 3.5, 6.5};
  std::vector<BinningValue> fCast = {Acts::binX, Acts::binY};
  auto fgAccess = GridAccessHelpers::castPosition<GridType>(gPosition, fCast);
  BOOST_CHECK_EQUAL(grid.atPosition(fgAccess), 300u);

  // Binned access
  Vector2 lbPosition = GridAccessHelpers::toLocal(grid, 4u, 9u);
  auto lbAccess =
      GridAccessHelpers::accessLocal<GridType>(lbPosition, fAccessor);
  BOOST_CHECK_EQUAL(grid.atPosition(lbAccess), 904u);
}

BOOST_AUTO_TEST_CASE(Grid3DAccess) {
  using EAxis = Acts::detail::Axis<Acts::detail::AxisType::Equidistant,
                                   Acts::detail::AxisBoundaryType::Bound>;
  using EGrid = Acts::Grid<std::size_t, EAxis, EAxis, EAxis>;

  auto xAxis = EAxis(0., 10., 10);
  auto yAxis = EAxis(0., 10., 10);
  auto zAxis = EAxis(0., 10., 10);

  EGrid grid({xAxis, yAxis, zAxis});

  BOOST_CHECK_THROW(GridAccessHelpers::accessLocal<EGrid>({1., 1.}, {0u, 1u});
                    , std::invalid_argument);

  BOOST_CHECK_THROW(GridAccessHelpers::toLocal(grid, 4u, 9u),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
