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

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/navigation/intersection/intersector_base.hpp"
#include "detray/navigation/intersection/ray_cylinder_intersector.hpp"
#include "detray/navigation/intersection/ray_cylinder_portal_intersector.hpp"
#include "detray/navigation/intersection/ray_line_intersector.hpp"
#include "detray/navigation/intersection/ray_plane_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_cylinder_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_cylinder_portal_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_line_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_plane_intersector.hpp"

namespace detray {

/// @brief Intersection implementation for detector surfaces using a ray
/// trajectory.
///
/// @note specialized into the concrete intersectors for the different local
/// geometries in the respective header files
template <typename frame_t, concepts::algebra algebra_t, bool resolve_pos>
struct ray_intersector_impl {};

template <typename shape_t, concepts::algebra algebra_t,
          bool resolve_pos = false>
using ray_intersector = intersector_base<
    ray_intersector_impl<typename shape_t::template local_frame_type<algebra_t>,
                         algebra_t, resolve_pos>>;

}  // namespace detray
