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
#include <cstdint>
#include <ostream>

namespace detray {

/// Shape of a detector volume.
///
/// cylinder: cylinder and two disc portals (shape: cylinder2D or cylinder3D).
/// cone: a cone portal and a disc portal (shape: missing).
/// rectangle: six rectangular portals that form a box (shape: cuboid3D).
/// trapezoid: six trapezoid portals (shape: cuboid3D).
/// cuboid: general cuboid form, excluding the previous ones (shape: cuboid3D).
enum class volume_id : std::uint_least8_t {
  e_cylinder = 0u,
  e_rectangle = 1u,
  e_trapezoid = 2u,
  e_cone = 3u,
  e_cuboid = 4u,
  e_unknown = 5u
};

/// surface type, resolved during navigation.
///
/// sensitive: can provide measurements and have material.
/// passive: no measurements, but can have material.
/// portal: boundary surface between two detector volumes, can have material.
enum class surface_id : std::uint_least8_t {
  e_portal = 0u,
  e_sensitive = 1u,
  e_passive = 2u,
  e_unknown = 3u,
  e_all = e_unknown
};

// Print the values of an enum by identifier
#define ENUM_PRINT(x) \
  case x:             \
    os << #x;         \
    break

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, volume_id vid) {
  switch (vid) {
    using enum volume_id;
    ENUM_PRINT(e_cylinder);
    ENUM_PRINT(e_rectangle);
    ENUM_PRINT(e_trapezoid);
    ENUM_PRINT(e_cone);
    ENUM_PRINT(e_cuboid);
    ENUM_PRINT(e_unknown);
  }
  return os;
}

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, surface_id sid) {
  switch (sid) {
    using enum surface_id;
    ENUM_PRINT(e_portal);
    ENUM_PRINT(e_sensitive);
    ENUM_PRINT(e_passive);
    case e_unknown:
      // e_all has same value (3u)
      os << "e_unknown/e_all";
      break;
  }
  return os;
}

#undef ENUM_PRINT
}  // namespace detray
