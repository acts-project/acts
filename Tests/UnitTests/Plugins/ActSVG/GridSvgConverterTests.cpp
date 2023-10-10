// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/ActSVG/GridSvgConverter.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"
#include "actsvg/display/grids.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

using namespace Acts;
using namespace Acts::detail;

namespace {
/// Helper method to turn a local bin into a string
///
/// @tparam local_bin_t the type for the local bin
/// @param lBin the local bin to printed
///
/// @return a string for screen output
template <typename local_bin_t, size_t DIM>
std::string localToString(const local_bin_t& lBin) {
  std::string lbString = "[";
  for (size_t ib = 0; ib < DIM; ++ib) {
    if (ib > 0u) {
      lbString += std::string(", ");
    }
    lbString += std::to_string(lBin[ib]);
  }
  lbString += std::string("]");
  return lbString;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(ActSvg)

BOOST_AUTO_TEST_CASE(BoundGridXY) {
  using GlobalBin = size_t;
  using LocalBin = std::array<size_t, 2u>;

  // x-y axis and grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisX(-200., 200, 4);
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisY(-200, 200, 6);
  Grid<std::tuple<GlobalBin, LocalBin>, decltype(axisX), decltype(axisY)>
      gridXY({axisX, axisY});

  Svg::GridConverter::Options cOptions;
  auto pGrid = Svg::GridConverter::convert(gridXY, {binX, binY}, cOptions);
  BOOST_CHECK(pGrid._type == actsvg::proto::grid::type::e_x_y);

  // Labelling the grid tiles
  auto edgesX = axisX.getBinEdges();
  auto edgesY = axisY.getBinEdges();

  std::vector<actsvg::svg::object> targets = {};
  for (auto [ix, x] : Acts::enumerate(edgesX)) {
    if (ix > 0u) {
      ActsScalar xp = 0.2 * edgesX[ix] + 0.8 * edgesX[ix - 1u];
      for (auto [iy, y] : Acts::enumerate(edgesY)) {
        if (iy > 0u) {
          ActsScalar yp = 0.8 * edgesY[iy] + 0.2 * edgesY[iy - 1u];
          decltype(gridXY)::point_t p = {xp, yp};
          // Get local and global index
          auto g = gridXY.globalBinFromPosition(p);
          auto l = gridXY.localBinsFromPosition(p);
          std::string gBin = std::string("g = ") + std::to_string(g);
          std::string lBin =
              std::string("l = ") +
              localToString<decltype(gridXY)::index_t, decltype(gridXY)::DIM>(
                  l);
          std::vector<std::string> glBin = {gBin, lBin};
          std::string gBinID = "g_" + std::to_string(g);
          targets.push_back(
              actsvg::draw::text(gBinID,
                                 {static_cast<actsvg::scalar>(xp),
                                  static_cast<actsvg::scalar>(yp)},
                                 glBin));
        }
      }
    }
  }
  pGrid._connections = targets;

  auto oGrid = actsvg::display::grid("BoundGridXY", pGrid);

  // Add some labelling
  actsvg::style::stroke axis_stroke{{{0, 0, 255}}, 3};
  actsvg::style::marker axis_marker{{"<"}, 4, {{{0, 0, 255}}}, axis_stroke};
  actsvg::style::font axis_font{{{0, 0, 255}}, "Andale Mono", 16};

  auto xAxis = actsvg::draw::arrow("x_axis", {0, 0}, {220, 0}, axis_stroke,
                                   actsvg::style::marker({""}), axis_marker);
  auto xLabel = actsvg::draw::text("x_label", {230, -4}, {"x"}, axis_font);
  auto yAxis = actsvg::draw::arrow("y_axis", {0, 0}, {0, 220}, axis_stroke,
                                   actsvg::style::marker({""}), axis_marker);
  auto yLabel = actsvg::draw::text("y_label", {-4, 232}, {"y"}, axis_font);

  std::vector<std::string> captionText = {
      "Binning schema for global and local bins: ",
      "- axis 0 : AxisBoundaryType::Bound, (-200., 200, 4), binX",
      "- axis 1 : AxisBoundaryType::Bound, (-200, 200, 6), binY"};

  auto caption = actsvg::draw::text("caption", {-180, -220}, captionText);

  oGrid.add_object(xAxis);
  oGrid.add_object(xLabel);
  oGrid.add_object(yAxis);
  oGrid.add_object(yLabel);
  oGrid.add_object(caption);

  Acts::Svg::toFile({oGrid}, oGrid._id + ".svg");
}

BOOST_AUTO_TEST_CASE(OpenGridXY) {
  using GlobalBin = size_t;
  using LocalBin = std::array<size_t, 2u>;

  // x-y axis and grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Open> axisX(-200., 200, 4);
  Axis<AxisType::Equidistant, AxisBoundaryType::Open> axisY(-200, 200, 6);
  Grid<std::tuple<GlobalBin, LocalBin>, decltype(axisX), decltype(axisY)>
      gridXY({axisX, axisY});

  Svg::GridConverter::Options cOptions;
  auto pGrid = Svg::GridConverter::convert(gridXY, {binX, binY}, cOptions);
  BOOST_CHECK(pGrid._type == actsvg::proto::grid::type::e_x_y);

  // Labelling the grid tiles
  auto edgesX = axisX.getBinEdges();
  auto edgesY = axisY.getBinEdges();

  std::vector<actsvg::svg::object> targets = {};
  size_t ig = 0;
  for (auto [ix, x] : Acts::enumerate(edgesX)) {
    if (ix > 0u) {
      ActsScalar xp = 0.2 * edgesX[ix] + 0.8 * edgesX[ix - 1u];
      for (auto [iy, y] : Acts::enumerate(edgesY)) {
        if (iy > 0u) {
          ActsScalar yp = 0.8 * edgesY[iy] + 0.2 * edgesY[iy - 1u];
          decltype(gridXY)::point_t p = {xp, yp};
          // Get local and global index
          auto g = gridXY.globalBinFromPosition(p);
          auto l = gridXY.localBinsFromPosition(p);
          std::string gBin = std::string("g = ") + std::to_string(g);
          std::string lBin =
              std::string("l = ") +
              localToString<decltype(gridXY)::index_t, decltype(gridXY)::DIM>(
                  l);
          std::vector<std::string> glBin = {gBin, lBin};
          std::string gBinID = "g_" + std::to_string(ig++);
          targets.push_back(
              actsvg::draw::text(gBinID,
                                 {static_cast<actsvg::scalar>(xp),
                                  static_cast<actsvg::scalar>(yp)},
                                 glBin));
        }
      }
    }
  }
  pGrid._connections = targets;

  // Add some labelling
  actsvg::style::stroke axis_stroke{{{0, 0, 255}}, 3};
  actsvg::style::marker axis_marker{{"<"}, 4, {{{0, 0, 255}}}, axis_stroke};
  actsvg::style::font axis_font{{{0, 0, 255}}, "Andale Mono", 16};

  auto xAxis = actsvg::draw::arrow("x_axis", {0, 0}, {220, 0}, axis_stroke,
                                   actsvg::style::marker({""}), axis_marker);
  auto xLabel = actsvg::draw::text("x_label", {230, -4}, {"x"}, axis_font);
  auto yAxis = actsvg::draw::arrow("y_axis", {0, 0}, {0, 220}, axis_stroke,
                                   actsvg::style::marker({""}), axis_marker);
  auto yLabel = actsvg::draw::text("y_label", {-4, 232}, {"y"}, axis_font);

  std::vector<std::string> captionText = {
      "Binning schema for global and local bins: ",
      "- axis 0 : AxisBoundaryType::Open, (-200., 200, 4), binX",
      "- axis 1 : AxisBoundaryType::Open, (-200, 200, 6), binY"};

  auto caption = actsvg::draw::text("caption", {-180, -220}, captionText);
  auto oGrid = actsvg::display::grid("OpenGridXY", pGrid);

  oGrid.add_object(xAxis);
  oGrid.add_object(xLabel);
  oGrid.add_object(yAxis);
  oGrid.add_object(yLabel);
  oGrid.add_object(caption);

  Acts::Svg::toFile({oGrid}, oGrid._id + ".svg");
}

BOOST_AUTO_TEST_CASE(ClosedCylinderGridZPhi) {
  using GlobalBin = size_t;
  using LocalBin = std::array<size_t, 2u>;

  // z-phi Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisZ(-200., 200., 3);
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(-M_PI, M_PI, 6);
  Grid<std::tuple<GlobalBin, LocalBin>, decltype(axisZ), decltype(axisPhi)>
      gridZPhi({axisZ, axisPhi});

  Svg::GridConverter::Options cOptions;
  auto pGrid = Svg::GridConverter::convert(gridZPhi, {binZ, binPhi}, cOptions);
  BOOST_CHECK(pGrid._type == actsvg::proto::grid::type::e_z_phi);

  pGrid._reference_r = 80.;

  // Labelling the grid tiles
  auto edgesZ = axisZ.getBinEdges();
  auto edgesPhi = axisPhi.getBinEdges();

  std::vector<actsvg::svg::object> targets = {};
  size_t ig = 0;
  for (auto [iz, z] : Acts::enumerate(edgesZ)) {
    if (iz > 0u) {
      ActsScalar zp = 0.2 * edgesZ[iz] + 0.8 * edgesZ[iz - 1u];
      for (auto [iphi, phi] : Acts::enumerate(edgesPhi)) {
        if (iphi > 0u) {
          ActsScalar phip = 0.8 * edgesPhi[iphi] + 0.2 * edgesPhi[iphi - 1u];
          decltype(gridZPhi)::point_t p = {zp, phip};
          // Get local and global index
          auto g = gridZPhi.globalBinFromPosition(p);
          auto l = gridZPhi.localBinsFromPosition(p);
          std::string gBin = std::string("g = ") + std::to_string(g);
          std::string lBin =
              std::string("l = ") + localToString<decltype(gridZPhi)::index_t,
                                                  decltype(gridZPhi)::DIM>(l);
          std::vector<std::string> glBin = {gBin, lBin};
          std::string gBinID = "g_" + std::to_string(ig++);
          targets.push_back(actsvg::draw::text(
              gBinID,
              {static_cast<actsvg::scalar>(zp),
               static_cast<actsvg::scalar>(pGrid._reference_r * phip)},
              glBin));
        }
      }
    }
  }
  pGrid._connections = targets;

  // Add some labelling
  actsvg::style::stroke axis_stroke{{{0, 0, 255}}, 3};
  actsvg::style::marker axis_marker{{"<"}, 4, {{{0, 0, 255}}}, axis_stroke};
  actsvg::style::font axis_font{{{0, 0, 255}}, "Andale Mono", 16};

  auto xAxis = actsvg::draw::arrow("x_axis", {0, 0}, {220, 0}, axis_stroke,
                                   actsvg::style::marker({""}), axis_marker);
  auto xLabel = actsvg::draw::text("x_label", {230, -4}, {"z"}, axis_font);
  auto yAxis = actsvg::draw::arrow("y_axis", {0, 0}, {0, 260}, axis_stroke,
                                   actsvg::style::marker({""}), axis_marker);
  auto yLabel = actsvg::draw::text("y_label", {-4, 270}, {"phi"}, axis_font);

  std::vector<std::string> captionText = {
      "Binning schema for global and local bins: ",
      "- axis 0 : AxisBoundaryType::Bound, (-200., 200, 3), binZ",
      "- axis 1 : AxisBoundaryType::Closed, (-PI, PI, 6), binPhi",
      "- draw reference radius set to 80"};

  auto caption = actsvg::draw::text("caption", {-180, -270}, captionText);
  auto oGrid = actsvg::display::grid("ClosedCylinderGridZPhi", pGrid);

  oGrid.add_object(xAxis);
  oGrid.add_object(xLabel);
  oGrid.add_object(yAxis);
  oGrid.add_object(yLabel);
  oGrid.add_object(caption);

  Acts::Svg::toFile({oGrid}, oGrid._id + ".svg");
}

BOOST_AUTO_TEST_CASE(ClosedDiscGridRPhi) {
  using GlobalBin = size_t;
  using LocalBin = std::array<size_t, 2u>;

  // r-phi Axes & Grid
  Axis<AxisType::Equidistant, AxisBoundaryType::Bound> axisR(100., 400., 3);
  Axis<AxisType::Equidistant, AxisBoundaryType::Closed> axisPhi(-M_PI, M_PI, 4);
  Grid<std::tuple<GlobalBin, LocalBin>, decltype(axisR), decltype(axisPhi)>
      gridRPhi({axisR, axisPhi});

  Svg::GridConverter::Options cOptions;
  auto pGrid = Svg::GridConverter::convert(gridRPhi, {binR, binPhi}, cOptions);
  BOOST_CHECK(pGrid._type == actsvg::proto::grid::type::e_r_phi);

  // Labelling the grid tiles
  auto edgesR = axisR.getBinEdges();
  auto edgesPhi = axisPhi.getBinEdges();

  std::vector<actsvg::svg::object> targets = {};
  size_t ig = 0;
  for (auto [ir, r] : Acts::enumerate(edgesR)) {
    if (ir > 0u) {
      ActsScalar rp = 0.5 * (edgesR[ir] + edgesR[ir - 1u]);
      for (auto [iphi, phi] : Acts::enumerate(edgesPhi)) {
        if (iphi > 0u) {
          ActsScalar phip = 0.5 * (edgesPhi[iphi] + edgesPhi[iphi - 1u]);
          decltype(gridRPhi)::point_t p = {rp, phip};
          // Get local and global index
          auto g = gridRPhi.globalBinFromPosition(p);
          auto l = gridRPhi.localBinsFromPosition(p);
          std::string gBin = std::string("g = ") + std::to_string(g);
          std::string lBin =
              std::string("l = ") + localToString<decltype(gridRPhi)::index_t,
                                                  decltype(gridRPhi)::DIM>(l);
          std::vector<std::string> glBin = {gBin, lBin};
          std::string gBinID = "g_" + std::to_string(ig++);
          targets.push_back(
              actsvg::draw::text(gBinID,
                                 {static_cast<actsvg::scalar>(rp * cos(phip)),
                                  static_cast<actsvg::scalar>(rp * sin(phip))},
                                 glBin));
        }
      }
    }
  }
  pGrid._connections = targets;

  // Add some labelling
  actsvg::style::stroke axis_stroke{{{0, 0, 255}}, 3};
  actsvg::style::marker axis_marker{{"<"}, 4, {{{0, 0, 255}}}, axis_stroke};
  actsvg::style::font axis_font{{{0, 0, 255}}, "Andale Mono", 16};

  auto rAxis = actsvg::draw::arrow("r_axis", {0, 0}, {420, 0}, axis_stroke,
                                   actsvg::style::marker({""}), axis_marker);
  auto rLabel = actsvg::draw::text("r_label", {440, -4}, {"r"}, axis_font);

  auto phiAxis = actsvg::draw::arc_measure(
      "phi_axis", 410., {410, 0.},
      {static_cast<actsvg::scalar>(410. * cos(0.25)),
       static_cast<actsvg::scalar>(410. * sin(0.25))},
      axis_stroke, actsvg::style::marker(), axis_marker);

  auto phiLabel =
      actsvg::draw::text("phi_label", {410, 60}, {"phi"}, axis_font);

  std::vector<std::string> captionText = {
      "Binning schema for global and local bins: ",
      "- axis 0 : AxisBoundaryType::Bound, (100., 400, 3), binR",
      "- axis 1 : AxisBoundaryType::Closed, (-PI, PI, 4), binPhi"};

  auto caption = actsvg::draw::text("caption", {-180, -420}, captionText);
  auto oGrid = actsvg::display::grid("ClosedDiscGridRPhi", pGrid);

  oGrid.add_object(rAxis);
  oGrid.add_object(rLabel);
  oGrid.add_object(phiAxis);
  oGrid.add_object(phiLabel);
  oGrid.add_object(caption);

  Acts::Svg::toFile({oGrid}, oGrid._id + ".svg");
}

BOOST_AUTO_TEST_SUITE_END()
