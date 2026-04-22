// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Projetc include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

namespace detray {

// some useful type declarations
using metadata_t = test::toy_metadata;
using detector_host_t = detector<metadata_t, host_container_types>;
using detector_device_t = detector<metadata_t, device_container_types>;

using det_volume_t = typename detector_host_t::volume_type;
using det_surface_t = typename detector_host_t::surface_type;
using transform_t = typename detector_host_t::transform3_type;
using mask_defs = typename detector_host_t::masks;

constexpr auto rectangle_id = mask_defs::id::e_rectangle2D;
constexpr auto disc_id = mask_defs::id::e_ring2D;
constexpr auto cylinder_id = mask_defs::id::e_concentric_cylinder2D;

using rectangle_t = types::get<mask_defs, rectangle_id>;
using disc_t = types::get<mask_defs, disc_id>;
using cylinder_t = types::get<mask_defs, cylinder_id>;

/// declaration of a test function for detector
void detector_test(typename detector_host_t::view_type det_data,
                   vecmem::data::vector_view<det_volume_t> volumes_data,
                   vecmem::data::vector_view<det_surface_t> surfaces_data,
                   vecmem::data::vector_view<transform_t> transforms_data,
                   vecmem::data::vector_view<rectangle_t> rectangles_data,
                   vecmem::data::vector_view<disc_t> discs_data,
                   vecmem::data::vector_view<cylinder_t> cylinders_data);

}  // namespace detray
