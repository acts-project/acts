// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <ostream>

namespace detray::axis {

/// axis bounds names.
///
/// open: adds additional over-/underflow bins.
/// closed: over-/underflow values are mapped into the bin range.
/// circular: over-/underflow values wrap around.
enum class bounds {
  e_open = 0,
  e_closed = 1,
  e_circular = 2,
};

/// axis coordinate names. Used to get a specific axes from an axes collection.
///
/// x, y, z: cartesian coordinate axes.
/// r, phi: polar coordinate axes.
/// rphi, cyl_z: 2D cylinder axes (3D cylinder uses r, phi, z).
enum class label {
  e_x = 0,
  e_y = 1,
  e_z = 2,
  e_r = 0,
  e_phi = 1,
  e_rphi = 0,
  e_cyl_z = 1,
};

/// axis binning type names.
///
/// regular: same sized bins along the axis.
/// irregular: every bin can have a different size along the axis.
enum class binning {
  e_regular = 0,
  e_irregular = 1,
};

// Print the values of an enum by identifier
#define ENUM_PRINT(x) \
  case x:             \
    os << #x;         \
    break

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, bounds b) {
  switch (b) {
    using enum bounds;
    ENUM_PRINT(e_open);
    ENUM_PRINT(e_closed);
    ENUM_PRINT(e_circular);
  }
  return os;
}

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, label l) {
  switch (l) {
    using enum label;
    case e_x:
      // e_r and e_rphi have same value (0)
      os << "e_x/e_r/e_rphi";
      break;
    case e_y:
      // e_phi and e_cyl_z have same value (1)
      os << "e_y/e_phi/e_cyl_z";
      break;
      ENUM_PRINT(e_z);
  }
  return os;
}

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, binning b) {
  switch (b) {
    using enum binning;
    ENUM_PRINT(e_regular);
    ENUM_PRINT(e_irregular);
  }
  return os;
}

#undef ENUM_PRINT

}  // namespace detray::axis
