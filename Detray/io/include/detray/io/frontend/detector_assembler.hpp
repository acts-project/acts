// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/io/backend/geometry_reader.hpp"
#include "detray/io/backend/homogeneous_material_reader.hpp"
#include "detray/io/backend/material_map_reader.hpp"
#include "detray/io/backend/surface_grid_reader.hpp"
#include "detray/io/frontend/definitions.hpp"
#include "detray/io/frontend/payloads.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace detray::io {

/// @brief Assemble a detray detector from its individual IO payloads.
///
/// This is the payload-driven counterpart of @c read_detector: instead of
/// reading the components from files, it consumes payloads that a client (e.g.
/// the ACTS geometry converter) has produced in memory and drives the detray
/// @c detector_builder to construct the detector. It is a single, extern-able
/// entry point for the otherwise heavy detector-building template tree.
///
/// @tparam detector_t the type of detector to be built
///
/// @param mr the memory resource used for the detector container allocations
/// @param geometry the geometry payload (always required)
/// @param homogeneous_material homogeneous material payload, or @c nullptr to
///        skip homogeneous material
/// @param material_grids material map payload (moved-from), or @c nullptr to
///        skip material maps
/// @param surface_grids surface grid payload, or @c nullptr to skip surface
///        grids
/// @param name detector name to set, or empty to leave the builder default
/// @param names filled with the volume/surface name map during the build
///
/// @returns the assembled detector
template <class detector_t>
detector_t assemble_detector(
    vecmem::memory_resource& mr, const detector_payload& geometry,
    const detector_homogeneous_material_payload* homogeneous_material,
    detector_grids_payload<surface_material_payload, material_id>*
        material_grids,
    const detector_grids_payload<std::size_t, accel_id>* surface_grids,
    std::string_view name, typename detector_t::name_map& names) {
  detector_builder<typename detector_t::metadata> det_builder{};

  geometry_reader::from_payload<detector_t>(det_builder, geometry);

  if (homogeneous_material != nullptr) {
    homogeneous_material_reader::from_payload<detector_t>(
        det_builder, *homogeneous_material);
  }

  if (material_grids != nullptr) {
    material_map_reader<std::integral_constant<std::size_t, 2>>::from_payload<
        detector_t>(det_builder, std::move(*material_grids));
  }

  if (surface_grids != nullptr) {
    surface_grid_reader<typename detector_t::surface_type,
                        std::integral_constant<std::size_t, 0>,
                        std::integral_constant<std::size_t, 2>>::
        template from_payload<detector_t>(det_builder, *surface_grids);
  }

  if (!name.empty()) {
    det_builder.set_name(std::string{name});
  }

  return det_builder.build(mr, names);
}

}  // namespace detray::io
