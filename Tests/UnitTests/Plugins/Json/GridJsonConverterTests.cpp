// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"

#include <array>
#include <fstream>
#include <memory>
#include <vector>

#include <nlohmann/json.hpp>

BOOST_AUTO_TEST_SUITE(GridJsonConversion)

BOOST_AUTO_TEST_CASE(Grid1DSingleEntry) {
  // Bound equidistant
  using EqBound = Acts::Experimental::detail::GridAxisGenerators::EqBound;

  EqBound eqBound{{0., 5.}, 5};
  // Create the grid with the provided axis generator
  using GridTypeEQB = typename EqBound::template grid_type<std::size_t>;
  GridTypeEQB eqBoundGrid(eqBound());

  eqBoundGrid.at(1u) = 1u;
  eqBoundGrid.at(2u) = 2u;
  eqBoundGrid.at(3u) = 3u;
  eqBoundGrid.at(4u) = 4u;
  eqBoundGrid.at(5u) = 5u;

  auto p1 = typename GridTypeEQB::point_t{0.5};
  BOOST_CHECK_EQUAL(eqBoundGrid.atPosition(p1), 1u);
  auto p2 = typename GridTypeEQB::point_t{1.5};
  BOOST_CHECK_EQUAL(eqBoundGrid.atPosition(p2), 2u);
  auto p3 = typename GridTypeEQB::point_t{2.5};
  BOOST_CHECK_EQUAL(eqBoundGrid.atPosition(p3), 3u);
  auto p4 = typename GridTypeEQB::point_t{3.5};
  BOOST_CHECK_EQUAL(eqBoundGrid.atPosition(p4), 4u);
  auto p5 = typename GridTypeEQB::point_t{4.5};
  BOOST_CHECK_EQUAL(eqBoundGrid.atPosition(p5), 5u);

  nlohmann::json eqBoundJson = Acts::GridJsonConverter::toJson(eqBoundGrid);

  auto eqBoundGridRead =
      Acts::GridJsonConverter::fromJson<EqBound, std::size_t>(eqBoundJson,
                                                              eqBound);

  BOOST_CHECK_EQUAL(eqBoundGridRead.at(1u), 1u);
  BOOST_CHECK_EQUAL(eqBoundGridRead.at(2u), 2u);
  BOOST_CHECK_EQUAL(eqBoundGridRead.at(3u), 3u);
  BOOST_CHECK_EQUAL(eqBoundGridRead.at(4u), 4u);
  BOOST_CHECK_EQUAL(eqBoundGridRead.at(5u), 5u);

  // Bound variable
  using VarBound = Acts::Experimental::detail::GridAxisGenerators::VarBound;

  VarBound varBound{{10., 11., 22., 333., 4444., 55555.}};
  // Create the grid with the provided axis generator
  using GridTypeEQV = typename VarBound::template grid_type<std::size_t>;
  GridTypeEQV varBoundGrid(varBound());

  varBoundGrid.at(1u) = 1u;
  varBoundGrid.at(2u) = 2u;
  varBoundGrid.at(3u) = 3u;
  varBoundGrid.at(4u) = 4u;
  varBoundGrid.at(5u) = 5u;

  nlohmann::json varBoundJson = Acts::GridJsonConverter::toJson(varBoundGrid);

  auto varBoundGridRead =
      Acts::GridJsonConverter::fromJson<VarBound, std::size_t>(varBoundJson,
                                                               varBound);

  BOOST_CHECK_EQUAL(varBoundGridRead.at(1u), 1u);
  BOOST_CHECK_EQUAL(varBoundGridRead.at(2u), 2u);
  BOOST_CHECK_EQUAL(varBoundGridRead.at(3u), 3u);
  BOOST_CHECK_EQUAL(varBoundGridRead.at(4u), 4u);
  BOOST_CHECK_EQUAL(varBoundGridRead.at(5u), 5u);

  // Closed equidistant
  using EqClosed = Acts::Experimental::detail::GridAxisGenerators::EqClosed;

  EqClosed eqClosed{{0., 5.}, 5};
  // Create the grid with the provided axis generator
  using GridTypeEQC = typename EqClosed::template grid_type<std::size_t>;
  GridTypeEQC eqClosedGrid(eqClosed());

  eqClosedGrid.at(1u) = 1u;
  eqClosedGrid.at(2u) = 2u;
  eqClosedGrid.at(3u) = 3u;
  eqClosedGrid.at(4u) = 4u;
  eqClosedGrid.at(5u) = 5u;

  nlohmann::json eqClosedJson = Acts::GridJsonConverter::toJson(eqClosedGrid);

  auto eqClosedGridRead =
      Acts::GridJsonConverter::fromJson<EqClosed, std::size_t>(eqClosedJson,
                                                               eqClosed);

  BOOST_CHECK_EQUAL(eqClosedGridRead.at(1u), 1u);
  BOOST_CHECK_EQUAL(eqClosedGridRead.at(2u), 2u);
  BOOST_CHECK_EQUAL(eqClosedGridRead.at(3u), 3u);
  BOOST_CHECK_EQUAL(eqClosedGridRead.at(4u), 4u);
  BOOST_CHECK_EQUAL(eqClosedGridRead.at(5u), 5u);
}

BOOST_AUTO_TEST_CASE(Grid1DArrayEntry) {
  // Bound equidistant
  using EqBound = Acts::Experimental::detail::GridAxisGenerators::EqBound;

  EqBound eqBound{{0., 5.}, 5};
  // Create the grid with the provided axis generator
  using GridTypeEQB =
      typename EqBound::template grid_type<std::array<std::size_t, 2u>>;
  GridTypeEQB eqBoundGrid(eqBound());

  eqBoundGrid.at(1u) = {1u, 1u};
  eqBoundGrid.at(2u) = {2u, 2u};
  eqBoundGrid.at(3u) = {3u, 3u};
  eqBoundGrid.at(4u) = {4u, 4u};
  eqBoundGrid.at(5u) = {5u, 5u};

  nlohmann::json eqBoundJson = Acts::GridJsonConverter::toJson(eqBoundGrid);

  auto eqBoundGridRead =
      Acts::GridJsonConverter::fromJson<EqBound, std::array<std::size_t, 2u>>(
          eqBoundJson, eqBound);

  BOOST_CHECK((eqBoundGridRead.at(1u) == std::array<std::size_t, 2u>{1u, 1u}));
  BOOST_CHECK((eqBoundGridRead.at(2u) == std::array<std::size_t, 2u>{2u, 2u}));
  BOOST_CHECK((eqBoundGridRead.at(3u) == std::array<std::size_t, 2u>{3u, 3u}));
  BOOST_CHECK((eqBoundGridRead.at(4u) == std::array<std::size_t, 2u>{4u, 4u}));
  BOOST_CHECK((eqBoundGridRead.at(5u) == std::array<std::size_t, 2u>{5u, 5u}));
}

BOOST_AUTO_TEST_CASE(Grid2DSingleEntryBound) {
  using EqBoundEqBound =
      Acts::Experimental::detail::GridAxisGenerators::EqBoundEqBound;

  EqBoundEqBound eqBound2{{0., 5.}, 5, {0., 2.}, 2};
  // Create the grid with the provided axis generator
  using GridTypeEQB2 = typename EqBoundEqBound::template grid_type<std::size_t>;
  GridTypeEQB2 eqBound2Grid(eqBound2());

  // Let's write in local coordinates
  using GridPoint = typename GridTypeEQB2::point_t;

  // First row access
  GridPoint p11{0.5, 0.5};
  GridPoint p12{1.5, 0.5};
  GridPoint p13{2.5, 0.5};
  GridPoint p14{3.5, 0.5};
  GridPoint p15{4.5, 0.5};
  eqBound2Grid.atPosition(p11) = 11u;
  eqBound2Grid.atPosition(p12) = 12u;
  eqBound2Grid.atPosition(p13) = 13u;
  eqBound2Grid.atPosition(p14) = 14u;
  eqBound2Grid.atPosition(p15) = 15u;

  // Second row access
  GridPoint p21{0.5, 1.5};
  GridPoint p22{1.5, 1.5};
  GridPoint p23{2.5, 1.5};
  GridPoint p24{3.5, 1.5};
  GridPoint p25{4.5, 1.5};
  eqBound2Grid.atPosition(p21) = 21u;
  eqBound2Grid.atPosition(p22) = 22u;
  eqBound2Grid.atPosition(p23) = 23u;
  eqBound2Grid.atPosition(p24) = 24u;
  eqBound2Grid.atPosition(p25) = 25u;

  nlohmann::json eqBound2Json = Acts::GridJsonConverter::toJson(eqBound2Grid);

  auto eqBound2JsonRead =
      Acts::GridJsonConverter::fromJson<EqBoundEqBound, std::size_t>(
          eqBound2Json, eqBound2);

  BOOST_CHECK_EQUAL(eqBound2JsonRead.atPosition(p11), 11u);
  BOOST_CHECK_EQUAL(eqBound2JsonRead.atPosition(p12), 12u);
  BOOST_CHECK_EQUAL(eqBound2JsonRead.atPosition(p13), 13u);
  BOOST_CHECK_EQUAL(eqBound2JsonRead.atPosition(p14), 14u);
  BOOST_CHECK_EQUAL(eqBound2JsonRead.atPosition(p15), 15u);
  BOOST_CHECK_EQUAL(eqBound2JsonRead.atPosition(p21), 21u);
  BOOST_CHECK_EQUAL(eqBound2JsonRead.atPosition(p22), 22u);
  BOOST_CHECK_EQUAL(eqBound2JsonRead.atPosition(p23), 23u);
  BOOST_CHECK_EQUAL(eqBound2JsonRead.atPosition(p24), 24u);
  BOOST_CHECK_EQUAL(eqBound2JsonRead.atPosition(p25), 25u);
}

BOOST_AUTO_TEST_CASE(Grid2DSingleEntryBoundClosed) {
  using EqBoundEqClosed =
      Acts::Experimental::detail::GridAxisGenerators::EqBoundEqClosed;

  EqBoundEqClosed eqBoundEqClosed{{-6., 6.}, 3, {-M_PI, M_PI}, 3};
  // Create the grid with the provided axis generator
  using GridTypeEQBEQC =
      typename EqBoundEqClosed::template grid_type<std::size_t>;
  GridTypeEQBEQC eqBoundEqClosedGrid(eqBoundEqClosed());

  // Let's write in local coordinates
  using GridPoint = typename GridTypeEQBEQC::point_t;

  // First row access
  GridPoint p11{-5, -2.};
  GridPoint p12{0., -2};
  GridPoint p13{5, -2};
  eqBoundEqClosedGrid.atPosition(p11) = 11u;
  eqBoundEqClosedGrid.atPosition(p12) = 12u;
  eqBoundEqClosedGrid.atPosition(p13) = 13u;

  // Middle row access
  GridPoint p21{-5., 0.};
  GridPoint p22{0., 0.};
  GridPoint p23{5., 0.};
  eqBoundEqClosedGrid.atPosition(p21) = 21u;
  eqBoundEqClosedGrid.atPosition(p22) = 22u;
  eqBoundEqClosedGrid.atPosition(p23) = 23u;

  // Last row access
  GridPoint p31{-5., 2.};
  GridPoint p32{0., 2.};
  GridPoint p33{5., 2.};
  eqBoundEqClosedGrid.atPosition(p31) = 31u;
  eqBoundEqClosedGrid.atPosition(p32) = 32u;
  eqBoundEqClosedGrid.atPosition(p33) = 33u;

  nlohmann::json eqBoundEqClosedJson =
      Acts::GridJsonConverter::toJson(eqBoundEqClosedGrid);

  auto eqBoundEqClosedJsonRead =
      Acts::GridJsonConverter::fromJson<EqBoundEqClosed, std::size_t>(
          eqBoundEqClosedJson, eqBoundEqClosed);

  BOOST_CHECK_EQUAL(eqBoundEqClosedJsonRead.atPosition(p11), 11u);
  BOOST_CHECK_EQUAL(eqBoundEqClosedJsonRead.atPosition(p12), 12u);
  BOOST_CHECK_EQUAL(eqBoundEqClosedJsonRead.atPosition(p13), 13u);

  BOOST_CHECK_EQUAL(eqBoundEqClosedJsonRead.atPosition(p21), 21u);
  BOOST_CHECK_EQUAL(eqBoundEqClosedJsonRead.atPosition(p22), 22u);
  BOOST_CHECK_EQUAL(eqBoundEqClosedJsonRead.atPosition(p23), 23u);

  BOOST_CHECK_EQUAL(eqBoundEqClosedJsonRead.atPosition(p31), 31u);
  BOOST_CHECK_EQUAL(eqBoundEqClosedJsonRead.atPosition(p32), 32u);
  BOOST_CHECK_EQUAL(eqBoundEqClosedJsonRead.atPosition(p33), 33u);
}

BOOST_AUTO_TEST_SUITE_END()
