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

// Project include(s)
#include "detray/material/concepts.hpp"
#include "detray/material/material_rod.hpp"
#include "detray/material/material_slab.hpp"
#include "detray/navigation/accelerators/concepts.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/type_registry.hpp"

// System include(s)
#include <type_traits>

namespace detray::detail {

/// Check for grid types in a data store
/// @{
template <typename>
struct contains_grids;

template <class ID, typename... Ts>
struct contains_grids<types::registry<ID, Ts...>> {
  static constexpr bool value{(concepts::grid<Ts> || ...)};
};

template <typename T>
inline constexpr bool contains_grids_v = contains_grids<T>::value;
/// @}

/// Check for surface grid types in a data store
/// @{
template <typename>
struct contains_surface_grids;

template <class ID, typename... Ts>
struct contains_surface_grids<types::registry<ID, Ts...>> {
  static constexpr bool value{(concepts::surface_grid<Ts> || ...)};
};

template <typename T>
inline constexpr bool contains_surface_grids_v =
    contains_surface_grids<T>::value;
/// @}

/// Check for the various types of material
/// @{

/// Contains slabs
/// @{
template <typename>
struct contains_material_slabs {};

template <class ID, typename... Ts>
struct contains_material_slabs<types::registry<ID, Ts...>> {
  static constexpr bool value{(concepts::material_slab<Ts> || ...)};
};

template <typename T>
inline constexpr bool contains_material_slabs_v =
    contains_material_slabs<T>::value;
/// @}

/// Contains rods
/// @{
template <typename>
struct contains_material_rods {};

template <class ID, typename... Ts>
struct contains_material_rods<types::registry<ID, Ts...>> {
  static constexpr bool value{(concepts::material_rod<Ts> || ...)};
};

template <typename T>
inline constexpr bool contains_material_rods_v =
    contains_material_rods<T>::value;
/// @}

/// Contains homogeneous material
/// @{
template <typename>
struct contains_homogeneous_material {};

template <class ID, typename... Ts>
struct contains_homogeneous_material<types::registry<ID, Ts...>> {
  static constexpr bool value{(concepts::homogeneous_material<Ts> || ...)};
};

template <typename T>
inline constexpr bool contains_homogeneous_material_v =
    contains_homogeneous_material<T>::value;
/// @}

/// Contains material maps
/// @{
template <typename>
struct contains_material_maps {};

template <class ID, typename... Ts>
struct contains_material_maps<types::registry<ID, Ts...>> {
  static constexpr bool value{(concepts::material_map<Ts> || ...)};
};

template <typename T>
inline constexpr bool contains_material_maps_v =
    contains_material_maps<T>::value;
/// @}

}  // namespace detray::detail
