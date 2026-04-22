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

namespace detray::navigation {

/// @enum NavigationDirection
/// The navigation direction is always with respect to a given track direction
enum class direction : std::int_least8_t { e_backward = -1, e_forward = 1 };

/// @enum Navigation status flags
enum class status : std::int_least8_t {
  e_abort = -3,          ///< error occurred, navigation will be aborted
  e_exit = -2,           ///< navigation finished/reached the end of geometry
  e_unknown = -1,        ///< unknown state/not initialized
  e_towards_object = 0,  ///< move towards next geometry object
  e_on_object = 1,       ///< reached a geometry object that is not a portal
  e_on_portal = 2,       ///< reached portal (material) surface
};

/// Navigation trust levels determine how the candidates cache is updated
enum class trust_level : std::uint_least8_t {
  e_no_trust = 0u,  ///< re-initialize the volume (i.e. run local navigation)
  e_fair = 1u,      ///< update the distance & order of the candidates
  e_high = 3u,      ///< update the dist. to the next candidate (current target)
  e_full = 4u       ///< don't update anything
};

// Print the values of an enum by identifier
#define ENUM_PRINT(x) \
  case x:             \
    os << #x;         \
    break

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, direction d) {
  switch (d) {
    using enum direction;
    ENUM_PRINT(e_backward);
    ENUM_PRINT(e_forward);
  }
  return os;
}

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, status st) {
  switch (st) {
    using enum status;
    ENUM_PRINT(e_abort);
    ENUM_PRINT(e_exit);
    ENUM_PRINT(e_unknown);
    ENUM_PRINT(e_towards_object);
    ENUM_PRINT(e_on_object);
    ENUM_PRINT(e_on_portal);
  }
  return os;
}

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, trust_level d) {
  switch (d) {
    using enum trust_level;
    ENUM_PRINT(e_no_trust);
    ENUM_PRINT(e_fair);
    ENUM_PRINT(e_high);
    ENUM_PRINT(e_full);
  }
  return os;
}

#undef ENUM_PRINT
}  // namespace detray::navigation
