// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/actor.hpp"

// System include(s)
#include <type_traits>

namespace detray {

/// Base class actor implementation
struct base_actor {
  /// Tag whether this is a composite type
  struct is_comp_actor : public std::false_type {};

  /// Defines the actors state. Hidden by actor implementations.
  struct state {};
};

namespace actor {

/// Dummy result type that signals that no result is present
struct empty_result {};

/// Result of the principal actor to be passed to the obsevers
struct result {
  actor::status status{actor::status::e_unknown};

  /// @returns a string stream that prints the transporter result details
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &os, const result &res) {
    os << "status: " << res.status << std::endl;
    return os;
  }
};

}  // namespace actor

}  // namespace detray
