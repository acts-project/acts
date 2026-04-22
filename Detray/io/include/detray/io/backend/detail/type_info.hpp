// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/io/frontend/definitions.hpp"
#include "detray/material/concepts.hpp"
#include "detray/navigation/accelerators/concepts.hpp"
#include "detray/utils/type_registry.hpp"

// System include(s)
#include <string_view>
#include <type_traits>

namespace detray::io::detail {

/// @brief empty type slot in registry (type unknown by given detector)
struct unknown_type {
  static constexpr std::string_view name = "unknown by detector";
};

/// Select mask shapes
/// @{
/// How to select shape types during construction for a given detector type
template <typename detector_t>
struct mask_shape_selector {
  using algebra_t = typename detector_t::algebra_type;
  using link_t = typename detector_t::surface_type::navigation_link;

  // Mask type in the detector corresponding to a given shape S
  template <typename S>
  using mask_t = detray::mask<S, algebra_t, link_t>;

  // The shape type is registered in the io shape_registry with a certain id
  template <typename S>
    requires concepts::shape<S, algebra_t>
  using type =
      std::conditional_t<types::contains<typename detector_t::masks, mask_t<S>>,
                         S, unknown_type>;
};

/// Type registry containing only the shapes required by the @tparam detector_t
template <typename detector_t>
using filtered_shape_registry =
    types::mapped_registry<io::shape_registry, mask_shape_selector<detector_t>>;

/// Infer the IO shape id from the shape type
template <typename shape_t>
  requires std::is_enum_v<typename shape_t::boundaries>
consteval io::shape_id get_id() {
  // Find the correct shape IO id;
  if constexpr (types::contains<io::shape_registry, shape_t>) {
    return types::id<io::shape_registry, shape_t>;
  } else {
    return io::shape_id::unknown;
  }
}
/// @}

/// Select material types
/// @{
/// Infer the IO material id from the material type - homogeneous material
template <detray::concepts::algebra algebra_t,
          detray::concepts::homogeneous_material material_t>
  requires std::same_as<dscalar<algebra_t>, typename material_t::scalar_type>
constexpr io::material_id get_id() {
  using mat_registry_t = io::material_registry<algebra_t>;

  // Find the correct material IO id;
  if constexpr (types::contains<mat_registry_t, material_t>) {
    return types::id<mat_registry_t, material_t>;
  } else {
    return io::material_id::unknown;
  }
}

/// Select the frame types corresponding to the local frames in the detector
/// material maps
/// @{
template <typename = void>
struct get_frame {
  using type = unknown_type;
};

template <concepts::material_map M>
struct get_frame<M> {
  using type = typename M::local_frame_type;
};

struct mat_map_frame_selector {
  template <typename M>
  using type = typename get_frame<M>::type;
};
/// @}

/// Type registry for the coord. frames required by the @tparam detector_t
template <typename detector_t>
using mat_frame_registry = types::mapped_registry<typename detector_t::material,
                                                  mat_map_frame_selector>;

/// How to select frame types during construction for a given detector type
template <typename detector_t>
struct mat_map_selector {
  template <typename F>
  using type =
      std::conditional_t<types::contains<mat_frame_registry<detector_t>, F>, F,
                         unknown_type>;
};

/// Type registry for the shapes required by the @tparam detector_t
template <typename detector_t>
using filtered_material_map_registry = types::mapped_registry<
    io::material_registry<typename detector_t::algebra_type>,
    mat_map_selector<detector_t>>;

/// Infer the IO material id from the material type - material maps
template <detray::concepts::material_map material_t>
constexpr io::material_id get_id() {
  using map_frame_t = typename material_t::local_frame_type;
  using algebra_t = typename map_frame_t::algebra_type;
  using mat_registry = io::material_registry<algebra_t>;

  // Find the correct material IO id;
  if constexpr (types::contains<mat_registry, map_frame_t>) {
    return types::id<mat_registry, map_frame_t>;
  } else {
    return io::material_id::unknown;
  }
}
/// @}

/// Select surface grid types
/// @{
/// Infer the grid id from its coordinate system
template <detray::concepts::surface_grid grid_t>
constexpr io::accel_id get_id() {
  using frame_t = typename grid_t::local_frame_type;
  using algebra_t = typename frame_t::algebra_type;
  using frame_registry_t = frame_registry<algebra_t>;

  // Find the correct grid shape IO id;
  if constexpr (types::contains<frame_registry_t, frame_t>) {
    return types::id<frame_registry_t, frame_t>;
  } else {
    return io::accel_id::unknown;
  }
}
/// @}

}  // namespace detray::io::detail
