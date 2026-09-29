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

/// Local navigation that a navigator requests from its caller in order to
/// complete an update. The caller decides when to run it.
enum class request : unsigned int {
  e_none = 0u,           ///< no local navigation needed
  e_init = 1u,           ///< initial local navigation
  e_re_init = 2u,        ///< re-initialize the current volume (trust lost)
  e_volume_switch = 3u,  ///< switch to the volume behind the current portal
  e_loose_re_init = 4u,  ///< re-initialize with loose tolerances (rescue)
};

/// Result of the cheap part of a navigation update
struct update_result {
  /// Local navigation that is needed to complete the update
  request next{request::e_none};
  /// Whether the navigation was re-initialized during the update itself
  bool is_init{false};
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
