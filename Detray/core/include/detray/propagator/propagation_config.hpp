// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/navigation/navigation_config.hpp"
#include "detray/propagator/stepping_config.hpp"

// System includes
#include <ostream>

namespace detray::propagation {

/// Configuration of the propagation
struct config {
  navigation::config navigation{};
  stepping::config stepping{};
  geometry_context context{};

  /// Print the propagation configuration
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& out, const config& cfg) {
    out << "Navigation\n"
        << "----------------------------\n"
        << cfg.navigation << "\nParameter Transport\n"
        << "----------------------------\n"
        << cfg.stepping << "\nGeometry Context\n"
        << "----------------------------\n"
        << "  No.                   : " << cfg.context.get() << "\n";

    return out;
  }
};

}  // namespace detray::propagation
