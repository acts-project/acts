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
#include "detray/definitions/units.hpp"

// System include(s)
#include <ostream>

namespace detray::intersection {

/// Intersector configuration
struct config {
  /// Tolerance on the mask 'is_inside' check:
  /// @{
  /// Minimal tolerance: ~ position uncertainty on surface
  float min_mask_tolerance{1e-5f * unit<float>::mm};
  /// Maximal tolerance: loose tolerance when still far away from surface
  float max_mask_tolerance{3.f * unit<float>::mm};
  /// Scale factor on the path used for the mask tolerance calculation
  float mask_tolerance_scalor{5e-2f};
  /// @}
  /// Maximal absolute path distance for a track to be considered 'on surface'
  float path_tolerance{1.f * unit<float>::um};
  /// How far behind the track position to look for candidates
  float overstep_tolerance{-1000.f * unit<float>::um};

  /// Print the intersector configuration
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& out, const config& cfg) {
    out << "  Min. mask tolerance   : "
        << cfg.min_mask_tolerance / detray::unit<float>::mm << " [mm]\n"
        << "  Max. mask tolerance   : "
        << cfg.max_mask_tolerance / detray::unit<float>::mm << " [mm]\n"
        << "  Mask tolerance scalor : " << cfg.mask_tolerance_scalor << "\n"
        << "  Path tolerance        : "
        << cfg.path_tolerance / detray::unit<float>::um << " [um]\n"
        << "  Overstep tolerance    : "
        << cfg.overstep_tolerance / detray::unit<float>::um << " [um]\n";

    return out;
  }
};

}  // namespace detray::intersection
