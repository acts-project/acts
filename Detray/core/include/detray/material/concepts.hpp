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

// Detray include(s)
#include "detray/material/material.hpp"
#include "detray/utils/grid/concepts.hpp"

namespace detray::concepts {

/// Material parameters
template <class M>
concept material_params = requires(const M m) {
  typename M::ratio;

  { m.X0() } -> std::same_as<typename M::scalar_type>;

  { m.L0() } -> std::same_as<typename M::scalar_type>;

  { m.Ar() } -> std::same_as<typename M::scalar_type>;

  { m.Z() } -> std::same_as<typename M::scalar_type>;

  { m.mass_density() } -> std::same_as<typename M::scalar_type>;

  { m.molar_density() } -> std::same_as<typename M::scalar_type>;

  { m.state() } -> std::same_as<detray::material_state>;

  { m.molar_electron_density() } -> std::same_as<typename M::scalar_type>;

  {
    m.density_effect_data()
  } -> std::same_as<
      const detray::detail::density_effect_data<typename M::scalar_type>&>;

  { m.mean_excitation_energy() } -> std::same_as<typename M::scalar_type>;
};

/// Material parameters with a thickness
template <class M>
concept material_slab = requires(const M slab) {
  typename M::material_type;

  requires concepts::material_params<typename M::material_type>;

  { slab.thickness() } -> std::same_as<typename M::scalar_type>;

  { slab.thickness_in_X0() } -> std::same_as<typename M::scalar_type>;

  { slab.thickness_in_L0() } -> std::same_as<typename M::scalar_type>;

  {
    slab.path_segment(typename M::scalar_type(), typename M::scalar_type())
  } -> std::same_as<typename M::scalar_type>;

  {
    slab.path_segment_in_X0(typename M::scalar_type(),
                            typename M::scalar_type())
  } -> std::same_as<typename M::scalar_type>;

  {
    slab.path_segment_in_L0(typename M::scalar_type(),
                            typename M::scalar_type())
  } -> std::same_as<typename M::scalar_type>;
};

/// Material parameters with a radius
template <class M>
concept material_rod = requires(const M rod) {
  typename M::material_type;

  requires concepts::material_params<typename M::material_type>;

  { rod.thickness() } -> std::same_as<typename M::scalar_type>;

  { rod.radius() } -> std::same_as<typename M::scalar_type>;

  {
    rod.path_segment(typename M::scalar_type(), typename M::scalar_type())
  } -> std::same_as<typename M::scalar_type>;

  {
    rod.path_segment_in_X0(typename M::scalar_type(), typename M::scalar_type())
  } -> std::same_as<typename M::scalar_type>;

  {
    rod.path_segment_in_L0(typename M::scalar_type(), typename M::scalar_type())
  } -> std::same_as<typename M::scalar_type>;
};

/// Homogeneous material
template <class M>
concept homogeneous_material =
    concepts::material_slab<M> || concepts::material_rod<M> ||
    concepts::material_params<M>;

/// Material map concept: Material grid
template <class M>
concept material_map =
    concepts::grid<M> && concepts::material_slab<typename M::value_type>;

/// Material that can used for surfaces
template <class M>
concept surface_material =
    (concepts::material_map<M> && (M::dim == 2)) ||
    concepts::material_slab<M> || concepts::material_rod<M>;

/// Material that can be used for volumes
template <class M>
concept volume_material = (concepts::material_map<M> && (M::dim == 3)) ||
                          concepts::material_params<M>;

}  // namespace detray::concepts
