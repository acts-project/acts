// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/math.hpp"
#include "detray/geometry/surface.hpp"

namespace detray {

/// Helper to get the incidence angle of a direction vector on a surface
///
/// @param ctx the geometric context
/// @param sf the geometric surface
/// @param dir normalized direction vector (e.g. direction of a track)
/// @param loc the local/bound position on the surface
///
/// @returns the cosine of the incidence angle given a local/bound position
template <typename detector_t, concepts::point point_t>
DETRAY_HOST_DEVICE constexpr dscalar<typename detector_t::algebra_type>
cos_angle(const typename detector_t::geometry_context &ctx,
          geometry::surface<detector_t> sf,
          const dvector3D<typename detector_t::algebra_type> &dir,
          const point_t &loc) {
  return math::fabs(vector::dot(dir, sf.normal(ctx, loc)));
}

}  // namespace detray
