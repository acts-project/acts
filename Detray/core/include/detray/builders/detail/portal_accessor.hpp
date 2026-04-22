// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/geometry/shapes/ring2D.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/geometry/tracking_volume.hpp"

// System include(s).
#include <vector>

namespace detray::detail {

/// @brief function that retrieves access to the portals of a cylindrical volume
///
/// The portals are returned as vectors in the order [inner, outer, lower,
/// upper]
template <typename detector_t>
auto get_cylinder_portals(const tracking_volume<detector_t> &vol) {
  using scalar_t = dscalar<typename detector_t::algebra_type>;

  std::vector<const typename detector_t::volume_type &> inner_pt{}, outer_pt{},
      lower_pt{}, upper_pr{};

  std::map<const typename detector_t::surface_type &, scalar_t> radii{0.f};
  std::map<const typename detector_t::surface_type &, scalar_t> z_pos{0.f};

  // Loop over all portals
  for (const auto &pt_desc : vol.portals()) {
    auto pt = geometry::surface{vol.detector(), pt_desc};
    const std::string name = pt.shape_name();

    if (name == "cylinder2D" || name == "concentric_cylinder2D") {
      radii[pt_desc] = pt.boundary(cylinder2D::e_r);
      z_pos[pt_desc] = pt.boundary(cylinder2D::e_lower_z);
      z_pos[pt_desc] = pt.boundary(cylinder2D::e_upper_z);
    } else {
      radii[pt_desc] = pt.boundary(ring2D::e_inner_r);
      radii[pt_desc] = pt.boundary(ring2D::e_outer_r);
    }
  }
}

}  // namespace detray::detail
