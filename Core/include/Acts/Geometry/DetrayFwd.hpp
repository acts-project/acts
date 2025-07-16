// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <variant>

namespace detray {

using dindex = unsigned int;

namespace io {

struct detector_payload;
struct transform_payload;
struct mask_payload;
struct surface_payload;
struct volume_payload;
struct material_slab_payload;
struct material_volume_payload;
struct detector_homogeneous_material_payload;

template <typename, typename>
struct grid_payload;
enum class material_id : unsigned int;

template <typename, typename>
struct detector_grids_payload;

enum class accel_id : unsigned int;

}  // namespace io
}  // namespace detray

namespace Acts {
using DetraySurfaceMaterial =
    std::variant<detray::io::grid_payload<detray::io::material_slab_payload,
                                          detray::io::material_id>,
                 detray::io::material_slab_payload>;

using DetraySurfaceGrid =
    detray::io::grid_payload<std::size_t, detray::io::accel_id>;

}  // namespace Acts
