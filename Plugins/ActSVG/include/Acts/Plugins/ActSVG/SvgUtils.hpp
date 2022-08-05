// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "actsvg/meta.hpp"

#include <array>
#include <fstream>
#include <string>
#include <vector>

namespace Acts {

namespace Svg {

struct Style {
  // Fill parameters
  std::array<int, 3> fillColor = {0, 0, 0};
  ActsScalar fillOpacity = 1.;

  // Highlight parameters
  std::array<int, 3> highlightColor = {0, 0, 0};
  std::vector<std::string> highlights = {};

  ActsScalar strokeWidth = 0.5;
  std::array<int, 3> strokeColor = {0, 0, 0};

  unsigned int nSegments = 72u;
};

/// Create a group
///
/// @param objects are the individual objects to be grouped
/// @param name is the name of the group
///
/// @return a signle svg object as a group
inline static actsvg::svg::object group(
    const std::vector<actsvg::svg::object>& objects, const std::string& name) {
  actsvg::svg::object gr;
  gr._tag = "g";
  gr._id = name;
  for (const auto& o : objects) {
    gr.add_object(o);
  }
  return gr;
}

/// Helper method to a measure
///
/// @param xStart the start position
/// @param yStart the start position
/// @param xEnd the start position
/// @param yEnd the start position
///
/// @return a single svg object as a measure
inline static actsvg::svg::object measure(ActsScalar xStart, ActsScalar yStart,
                                          ActsScalar xEnd, ActsScalar yEnd,
                                          const std::string& variable = "",
                                          ActsScalar value = 0.,
                                          const std::string& unit = "") {
  std::string mlabel = "";
  if (not variable.empty()) {
    mlabel = variable + " = ";
  }
  if (value != 0.) {
    mlabel += actsvg::utils::to_string(static_cast<actsvg::scalar>(value));
  }
  if (not unit.empty()) {
    mlabel += " ";
    mlabel += unit;
  }
  return actsvg::draw::measure(
      "measure",
      {static_cast<actsvg::scalar>(xStart),
       static_cast<actsvg::scalar>(yStart)},
      {static_cast<actsvg::scalar>(xEnd), static_cast<actsvg::scalar>(yEnd)},
      actsvg::style::stroke(), actsvg::style::marker({"o"}),
      actsvg::style::marker({"|<<"}), actsvg::style::font(), mlabel);
}

// Helper method to draw axes
///
/// @param xMin the minimum x value
/// @param xMax the maximum x value
/// @param yMin the minimum y value
/// @param yMax the maximum y value
///
/// @retrun an svg object
inline static actsvg::svg::object axesXY(ActsScalar xMin, ActsScalar xMax,
                                         ActsScalar yMin, ActsScalar yMax) {
  return actsvg::draw::x_y_axes(
      "x_y_axis",
      {static_cast<actsvg::scalar>(xMin), static_cast<actsvg::scalar>(xMax)},
      {static_cast<actsvg::scalar>(yMin), static_cast<actsvg::scalar>(yMax)});
}

/// Helper method to write to file
///
/// @param objects to be written out
/// @param fileName the file name is to be given
///
inline static void toFile(const std::vector<actsvg::svg::object>& objects,
                          const std::string& fileName) {
  actsvg::svg::file foutFile;

  for (const auto& o : objects) {
    foutFile.add_object(o);
  }

  std::ofstream fout;
  fout.open(fileName);
  fout << foutFile;
  fout.close();
}

}  // namespace Svg

}  // namespace Acts