// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/type_traits.hpp"

namespace detray::concepts {

/// Check for the the presence of any type of grids in a detector definition
template <class detector_t>
concept has_grids = detail::contains_grids_v<typename detector_t::accel> ||
                    detail::contains_grids_v<typename detector_t::material>;

/// Check for the the presence of surface grids in a detector definition
template <class detector_t>
concept has_surface_grids =
    detail::contains_surface_grids_v<typename detector_t::accel>;

/// Check for the the presence of material slabs in a detector definition
template <class detector_t>
concept has_material_slabs =
    detail::contains_material_slabs_v<typename detector_t::material>;

/// Check for the the presence of material rods in a detector definition
template <class detector_t>
concept has_material_rods =
    detail::contains_material_rods_v<typename detector_t::material>;

/// Check for the the presence of homogeneous material types in a detector
/// definition
template <class detector_t>
concept has_homogeneous_material =
    detail::contains_homogeneous_material_v<typename detector_t::material>;

/// Check for the the presence of material maps in a detector definition
template <class detector_t>
concept has_material_maps =
    detail::contains_material_maps_v<typename detector_t::material>;

}  // namespace detray::concepts
