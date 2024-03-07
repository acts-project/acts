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
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

using namespace Acts;

namespace bd = boost::unit_test::data;

using LocalAccess = GridAccessHelpers::LocalAccess;

BOOST_AUTO_TEST_SUITE(GridAccessHelpersTests)

BOOST_AUTO_TEST_CASE(LocalAccess_test) {
  // Pick out the x coordinate
  LocalAccess lAccess{0u};
  Vector2 position{3., 4.};

  BOOST_CHECK_EQUAL(lAccess.toGridLocal(position), 3.);

  // Assuma a cylindrical surface with a phi - z grid on it
  ActsScalar radius = 100;
  ActsScalar phiValue = 0.25;
  ActsScalar zValue = 55.;

  LocalAccess rphiToPhi{0u, 0., 1. / radius};
  LocalAccess zToZ{1u};
  std::vector<LocalAccess> lAccessors = {rphiToPhi, zToZ};

  using GridType =
      Acts::GridAxisGenerators::EqClosedEqBound::grid_type<std::size_t>;

  Vector2 rphiPosition{radius * phiValue, zValue};
  auto rphiAccess =
      GridAccessHelpers::castLocal<GridType>(rphiPosition, lAccessors);

  BOOST_CHECK_EQUAL(rphiAccess.size(), 2u);
  BOOST_CHECK_EQUAL(rphiAccess[0u], phiValue);
  BOOST_CHECK_EQUAL(rphiAccess[1u], zValue);

  // Check the reverse
  ActsScalar zFrame = zToZ.toFrameLocal(zValue);
  BOOST_CHECK_EQUAL(zFrame, zValue);

  ActsScalar loc0Frame = rphiToPhi.toFrameLocal(phiValue);
  CHECK_CLOSE_ABS(loc0Frame, phiValue * radius, 1e-6);

  // Assume the same cylidner, but someone maps a z - phi grid onto it
  std::vector<LocalAccess> lAccessorsRev = {zToZ, rphiToPhi};
  auto rphiAccessRev =
      GridAccessHelpers::castLocal<GridType>(rphiPosition, lAccessorsRev);
  BOOST_CHECK_EQUAL(rphiAccessRev.size(), 2u);
  BOOST_CHECK_EQUAL(rphiAccessRev[0u], zValue);
  BOOST_CHECK_EQUAL(rphiAccessRev[1u], phiValue);
}

BOOST_AUTO_TEST_CASE(Grid1DAccess) {
  Acts::GridAxisGenerators::EqBound eqBound{{0., 10.}, 10};
  using GridType = Acts::GridAxisGenerators::EqBound::grid_type<std::size_t>;
  using PointType = GridType::point_t;
  auto grid = GridType(eqBound());

  for (std::size_t i = 0; i < 10; ++i) {
    grid.atPosition(PointType{i + 0.5}) = i;
  }

  // Local access
  std::vector<LocalAccess> fAccessor = {LocalAccess{0u}};
  std::vector<LocalAccess> sAccessor = {LocalAccess{1u}};

  Vector2 lPosition{3.5, 6.5};
  auto flAccess = GridAccessHelpers::castLocal<GridType>(lPosition, fAccessor);
  auto slAccess = GridAccessHelpers::castLocal<GridType>(lPosition, sAccessor);

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
  std::vector<LocalAccess> fAccessor = {LocalAccess{0u}, LocalAccess{1u}};
  Vector2 lPosition{3.5, 6.5};
  auto flAccess = GridAccessHelpers::castLocal<GridType>(lPosition, fAccessor);
  BOOST_CHECK_EQUAL(grid.atPosition(flAccess), 603u);

  // Global access
  Vector3 gPosition{0.5, 3.5, 6.5};
  std::vector<BinningValue> fCast = {Acts::binX, Acts::binY};
  auto fgAccess = GridAccessHelpers::castPosition<GridType>(gPosition, fCast);
  BOOST_CHECK_EQUAL(grid.atPosition(fgAccess), 300u);
}

BOOST_AUTO_TEST_SUITE_END()
