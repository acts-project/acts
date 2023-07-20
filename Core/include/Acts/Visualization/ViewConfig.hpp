// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <string>

namespace Acts {

/// The color typedef. It's an array of three numbers [0, 255]
/// representing the RGB color values.
///
using ColorRGB = std::array<int, 3>;

/// @brief Struct to concentrate all visualization configurations
/// in order to harmonize visualization interfaces
struct ViewConfig {
  /// Constructor to switch visibility off
  ViewConfig(bool vis = true) : visible(vis) {}

  /// Constructor for color settings only
  ViewConfig(const ColorRGB& rgb) : color(rgb) {}

  /// Visible flag
  bool visible = true;
  /// The RGB color for this object
  ColorRGB color = {250, 0, 0};
  /// Out of plane drawing parameter for objects
  double offset = 0.1;
  /// The visual line thickness for this object
  double lineThickness = 0.15;
  /// The visual surface thickness for this object
  double surfaceThickness = 0.15;
  /// The number of segments to approximate full 2pi
  unsigned int nSegments = 72;
  /// Whether to triangulate or not
  bool triangulate = false;
  /// Write name - non-empty string indicates writing
  std::string outputName = "";
};

}  // namespace Acts
