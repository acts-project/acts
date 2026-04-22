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
#include "detray/core/concepts.hpp"
#include "detray/io/backend/geometry_writer.hpp"
#include "detray/io/backend/homogeneous_material_writer.hpp"
#include "detray/io/backend/material_map_writer.hpp"
#include "detray/io/backend/surface_grid_writer.hpp"
#include "detray/io/frontend/detail/detector_components_writer.hpp"
#include "detray/io/frontend/detector_writer_config.hpp"
#include "detray/io/json/json_converter.hpp"

namespace detray::io {

struct detector_writer_config;

namespace detail {

/// Infer the writers that are needed from the detector type @tparam detector_t
template <class detector_t>
void add_json_writers(detector_components_writer<detector_t>& writers,
                      const detray::io::detector_writer_config& cfg) {
  // Always needed
  using json_geometry_writer = json_converter<detector_t, geometry_writer>;

  writers.template add<json_geometry_writer>();

  // Find other writers, depending on the detector type
  if (cfg.write_material()) {
    // Simple material
    if constexpr (detray::concepts::has_homogeneous_material<detector_t>) {
      using json_homogeneous_material_writer =
          json_converter<detector_t, homogeneous_material_writer>;

      writers.template add<json_homogeneous_material_writer>();
    }
    // Material maps
    if constexpr (detray::concepts::has_material_maps<detector_t>) {
      using json_material_map_writer =
          json_converter<detector_t, material_map_writer>;

      writers.template add<json_material_map_writer>();
    }
  }
  // Navigation acceleration structures
  if constexpr (detray::concepts::has_surface_grids<detector_t>) {
    using json_surface_grid_writer =
        json_converter<detector_t, surface_grid_writer>;

    if (cfg.write_grids()) {
      writers.template add<json_surface_grid_writer>();
    }
  }
}

}  // namespace detail

}  // namespace detray::io
