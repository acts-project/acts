// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Actsvg include(s)
#include "actsvg/core.hpp"

// System include(s)
#include <vector>

namespace detray::svgtools::meta::proto {

/// @brief A proto eta line
struct eta_lines {
  actsvg::scalar _r{200.f};
  actsvg::scalar _z{800.f};

  // Main eta lines
  std::vector<actsvg::scalar> _values_main = {1.f, 2.f, 3.f, 4.f};
  // Intermediate eta lines
  std::vector<actsvg::scalar> _values_half = {0.5f, 1.5f, 2.5f, 3.5f};

  actsvg::style::stroke _stroke_main;
  actsvg::style::stroke _stroke_half;
  actsvg::style::font _label_font;
  bool _show_label = true;
};
}  // namespace detray::svgtools::meta::proto
