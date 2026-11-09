// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <vector>

namespace detray::svgtools::meta::proto {

/// @brief A proto eta line
struct eta_lines {
  /// Maximum radial extent of the eta-line overlay.
  actsvg::scalar _r{200.f};
  /// Maximum longitudinal extent of the eta-line overlay.
  actsvg::scalar _z{800.f};

  // Main eta lines
  /// Absolute eta values for the primary lines.
  std::vector<actsvg::scalar> _values_main = {1.f, 2.f, 3.f, 4.f};
  // Intermediate eta lines
  /// Absolute eta values for intermediate lines.
  std::vector<actsvg::scalar> _values_half = {0.5f, 1.5f, 2.5f, 3.5f};

  /// Stroke style for primary eta lines.
  actsvg::style::stroke _stroke_main;
  /// Stroke style for intermediate eta lines.
  actsvg::style::stroke _stroke_half;
  /// Font used for eta labels.
  actsvg::style::font _label_font;
  /// Whether eta values are drawn beside the lines.
  bool _show_label = true;
};
}  // namespace detray::svgtools::meta::proto
