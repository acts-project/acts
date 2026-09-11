// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/GridAxisGenerators.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"

#include <array>
#include <fstream>
#include <memory>
#include <numbers>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(Grid1DSingleEntry) {
  // Bound equidistant
  using EqBound = GridAxisGenerators::EqBound;

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

  nlohmann::json eqBoundJson = GridJsonConverter::toJson(eqBoundGrid);

  auto eqBoundGridRead =
      GridJsonConverter::fromJson<EqBound, std::size_t>(eqBoundJson, eqBound);

  BOOST_CHECK_EQUAL(eqBoundGridRead.at(1u), 1u);
  BOOST_CHECK_EQUAL(eqBoundGridRead.at(2u), 2u);
  BOOST_CHECK_EQUAL(eqBoundGridRead.at(3u), 3u);
  BOOST_CHECK_EQUAL(eqBoundGridRead.at(4u), 4u);
  BOOST_CHECK_EQUAL(eqBoundGridRead.at(5u), 5u);

  // Bound variable
  using VarBound = GridAxisGenerators::VarBound;

  VarBound varBound{{10., 11., 22., 333., 4444., 55555.}};
  // Create the grid with the provided axis generator
  using GridTypeEQV = typename VarBound::template grid_type<std::size_t>;
  GridTypeEQV varBoundGrid(varBound());

  varBoundGrid.at(1u) = 1u;
  varBoundGrid.at(2u) = 2u;
  varBoundGrid.at(3u) = 3u;
  varBoundGrid.at(4u) = 4u;
  varBoundGrid.at(5u) = 5u;

  nlohmann::json varBoundJson = GridJsonConverter::toJson(varBoundGrid);

  auto varBoundGridRead = GridJsonConverter::fromJson<VarBound, std::size_t>(
      varBoundJson, varBound);

  BOOST_CHECK_EQUAL(varBoundGridRead.at(1u), 1u);
  BOOST_CHECK_EQUAL(varBoundGridRead.at(2u), 2u);
  BOOST_CHECK_EQUAL(varBoundGridRead.at(3u), 3u);
  BOOST_CHECK_EQUAL(varBoundGridRead.at(4u), 4u);
  BOOST_CHECK_EQUAL(varBoundGridRead.at(5u), 5u);

  // Closed equidistant
  using EqClosed = GridAxisGenerators::EqClosed;

  EqClosed eqClosed{{0., 5.}, 5};
  // Create the grid with the provided axis generator
  using GridTypeEQC = typename EqClosed::template grid_type<std::size_t>;
  GridTypeEQC eqClosedGrid(eqClosed());

  eqClosedGrid.at(1u) = 1u;
  eqClosedGrid.at(2u) = 2u;
  eqClosedGrid.at(3u) = 3u;
  eqClosedGrid.at(4u) = 4u;
  eqClosedGrid.at(5u) = 5u;

  nlohmann::json eqClosedJson = GridJsonConverter::toJson(eqClosedGrid);

  auto eqClosedGridRead = GridJsonConverter::fromJson<EqClosed, std::size_t>(
      eqClosedJson, eqClosed);

  BOOST_CHECK_EQUAL(eqClosedGridRead.at(1u), 1u);
  BOOST_CHECK_EQUAL(eqClosedGridRead.at(2u), 2u);
  BOOST_CHECK_EQUAL(eqClosedGridRead.at(3u), 3u);
  BOOST_CHECK_EQUAL(eqClosedGridRead.at(4u), 4u);
  BOOST_CHECK_EQUAL(eqClosedGridRead.at(5u), 5u);
}

BOOST_AUTO_TEST_CASE(Grid1DArrayEntry) {
  // Bound equidistant
  using EqBound = GridAxisGenerators::EqBound;

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

  nlohmann::json eqBoundJson = GridJsonConverter::toJson(eqBoundGrid);

  auto eqBoundGridRead =
      GridJsonConverter::fromJson<EqBound, std::array<std::size_t, 2u>>(
          eqBoundJson, eqBound);

  BOOST_CHECK((eqBoundGridRead.at(1u) == std::array<std::size_t, 2u>{1u, 1u}));
  BOOST_CHECK((eqBoundGridRead.at(2u) == std::array<std::size_t, 2u>{2u, 2u}));
  BOOST_CHECK((eqBoundGridRead.at(3u) == std::array<std::size_t, 2u>{3u, 3u}));
  BOOST_CHECK((eqBoundGridRead.at(4u) == std::array<std::size_t, 2u>{4u, 4u}));
  BOOST_CHECK((eqBoundGridRead.at(5u) == std::array<std::size_t, 2u>{5u, 5u}));
}

BOOST_AUTO_TEST_CASE(Grid2DSingleEntryBound) {
  using EqBoundEqBound = GridAxisGenerators::EqBoundEqBound;

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

  nlohmann::json eqBound2Json = GridJsonConverter::toJson(eqBound2Grid);

  auto eqBound2JsonRead =
      GridJsonConverter::fromJson<EqBoundEqBound, std::size_t>(eqBound2Json,
                                                               eqBound2);

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
  using EqBoundEqClosed = GridAxisGenerators::EqBoundEqClosed;

  EqBoundEqClosed eqBoundEqClosed{
      {-6., 6.}, 3, {-std::numbers::pi, std::numbers::pi}, 3};
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
      GridJsonConverter::toJson(eqBoundEqClosedGrid);

  auto eqBoundEqClosedJsonRead =
      GridJsonConverter::fromJson<EqBoundEqClosed, std::size_t>(
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

BOOST_AUTO_TEST_CASE(GridAnyToJson1D) {
  using EqBound = GridAxisGenerators::EqBound;

  EqBound eqBound{{0., 5.}, 5};
  using GridTypeEQB = typename EqBound::template grid_type<std::size_t>;
  GridTypeEQB eqBoundGrid(eqBound());

  eqBoundGrid.at(1u) = 1u;
  eqBoundGrid.at(2u) = 2u;
  eqBoundGrid.at(3u) = 3u;
  eqBoundGrid.at(4u) = 4u;
  eqBoundGrid.at(5u) = 5u;

  nlohmann::json jExpected = GridJsonConverter::toJson(eqBoundGrid);
  nlohmann::json jAny = GridJsonConverter::toJsonAny<std::size_t>(
      eqBoundGrid, AnyGridConstView<std::size_t>(eqBoundGrid));

  BOOST_CHECK_EQUAL(jAny, jExpected);
}

BOOST_AUTO_TEST_CASE(GridAnyToJson2D) {
  using EqBoundEqClosed = GridAxisGenerators::EqBoundEqClosed;

  EqBoundEqClosed eqBoundEqClosed{
      {-6., 6.}, 3, {-std::numbers::pi, std::numbers::pi}, 3};
  using GridTypeEQBEQC =
      typename EqBoundEqClosed::template grid_type<std::size_t>;
  GridTypeEQBEQC eqBoundEqClosedGrid(eqBoundEqClosed());

  using GridPoint = typename GridTypeEQBEQC::point_t;
  eqBoundEqClosedGrid.atPosition(GridPoint{-5., -2.}) = 11u;
  eqBoundEqClosedGrid.atPosition(GridPoint{0., -2.}) = 12u;
  eqBoundEqClosedGrid.atPosition(GridPoint{5., -2.}) = 13u;
  eqBoundEqClosedGrid.atPosition(GridPoint{-5., 0.}) = 21u;
  eqBoundEqClosedGrid.atPosition(GridPoint{0., 0.}) = 22u;
  eqBoundEqClosedGrid.atPosition(GridPoint{5., 0.}) = 23u;
  eqBoundEqClosedGrid.atPosition(GridPoint{-5., 2.}) = 31u;
  eqBoundEqClosedGrid.atPosition(GridPoint{0., 2.}) = 32u;
  eqBoundEqClosedGrid.atPosition(GridPoint{5., 2.}) = 33u;

  nlohmann::json jExpected = GridJsonConverter::toJson(eqBoundEqClosedGrid);
  nlohmann::json jAny = GridJsonConverter::toJsonAny<std::size_t>(
      eqBoundEqClosedGrid, AnyGridConstView<std::size_t>(eqBoundEqClosedGrid));

  BOOST_CHECK_EQUAL(jAny, jExpected);
}

BOOST_AUTO_TEST_CASE(AxisJsonConverterEquidistantBound) {
  auto axis = IAxis::createEquidistant(AxisBoundaryType::Bound, -5., 5., 10);
  nlohmann::json j = AxisJsonConverter::toJson(*axis);
  auto read = AxisJsonConverter::fromJson(j);
  BOOST_REQUIRE(read != nullptr);
  BOOST_CHECK_EQUAL(read->getBoundaryType(), AxisBoundaryType::Bound);
  BOOST_CHECK_EQUAL(read->getNBins(), 10u);
  BOOST_CHECK_EQUAL(read->getMin(), -5.);
  BOOST_CHECK_EQUAL(read->getMax(), 5.);
  BOOST_CHECK(read->isEquidistant());
}

BOOST_AUTO_TEST_CASE(AxisJsonConverterEquidistantClosed) {
  auto axis = IAxis::createEquidistant(AxisBoundaryType::Closed,
                                       -std::numbers::pi, std::numbers::pi, 36);
  nlohmann::json j = AxisJsonConverter::toJson(*axis);
  auto read = AxisJsonConverter::fromJson(j);
  BOOST_REQUIRE(read != nullptr);
  BOOST_CHECK_EQUAL(read->getBoundaryType(), AxisBoundaryType::Closed);
  BOOST_CHECK_EQUAL(read->getNBins(), 36u);
  BOOST_CHECK(read->isEquidistant());
}

BOOST_AUTO_TEST_CASE(AxisJsonConverterVariableBound) {
  // Exercises the "boundaries" key — the fix changed from reading "edges"
  std::vector<double> edges = {0., 10., 25., 60., 130.};
  auto axis = IAxis::createVariable(AxisBoundaryType::Bound, edges);
  nlohmann::json j = AxisJsonConverter::toJson(*axis);
  BOOST_CHECK_EQUAL(j.at("type").get<AxisType>(), AxisType::Variable);
  BOOST_CHECK(j.contains("boundaries"));
  BOOST_CHECK(!j.contains("edges"));
  auto read = AxisJsonConverter::fromJson(j);
  BOOST_REQUIRE(read != nullptr);
  BOOST_CHECK_EQUAL(read->getBoundaryType(), AxisBoundaryType::Bound);
  BOOST_CHECK_EQUAL(read->getNBins(), edges.size() - 1);
  BOOST_CHECK(!read->isEquidistant());
  BOOST_CHECK(read->getBinEdges() == edges);
}

BOOST_AUTO_TEST_CASE(AxisJsonConverterDirection) {
  // The axis direction is round-tripped when present ...
  auto axis = IAxis::createEquidistant(AxisBoundaryType::Bound, -5., 5., 10,
                                       AxisDirection::AxisZ);
  nlohmann::json j = AxisJsonConverter::toJson(*axis);
  BOOST_CHECK(j.contains("direction"));
  auto read = AxisJsonConverter::fromJson(j);
  BOOST_REQUIRE(read != nullptr);
  BOOST_CHECK(read->getDirection() == AxisDirection::AxisZ);

  // ... and stays absent when not set
  auto axisNoDir =
      IAxis::createEquidistant(AxisBoundaryType::Bound, -5., 5., 10);
  nlohmann::json jNoDir = AxisJsonConverter::toJson(*axisNoDir);
  BOOST_CHECK(!jNoDir.contains("direction"));
  auto readNoDir = AxisJsonConverter::fromJson(jNoDir);
  BOOST_REQUIRE(readNoDir != nullptr);
  BOOST_CHECK(!readNoDir->getDirection().has_value());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
