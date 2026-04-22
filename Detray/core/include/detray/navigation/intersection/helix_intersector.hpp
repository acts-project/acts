// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/navigation/intersection/helix_cylinder_intersector.hpp"
#include "detray/navigation/intersection/helix_line_intersector.hpp"
#include "detray/navigation/intersection/helix_plane_intersector.hpp"
#include "detray/navigation/intersection/intersector_base.hpp"

namespace detray {

/// @brief Intersection implementation for detector surfaces using a helix
/// trajectory.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
///
/// @note specialized into @c helix_plane_intersector, @c helix_line_intersector
/// and @c helix_cylinder_intersector
template <typename frame_t, concepts::algebra algebra_t>
struct helix_intersector_impl {};

template <typename shape_t, concepts::algebra algebra_t, bool = true>
using helix_intersector = intersector_base<helix_intersector_impl<
    typename shape_t::template local_frame_type<algebra_t>, algebra_t>>;

}  // namespace detray
