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
#include "detray/builders/detail/radius_getter.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/utils/grid/concepts.hpp"

// Plugin include(s)
#include "detray/plugins/svgtools/conversion/grid.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"

// Actsvg include(s)
#include "actsvg/proto/surface.hpp"

// System include(s)
#include <algorithm>
#include <iterator>

namespace detray::svgtools::conversion {

/// @brief Converts a the material grid of a detray surface to an actsvg proto
///        grid.
///
/// @param detector the detector
/// @param index the index of the surface
/// @param view the view
/// @param style the style settings
///
/// @returns a proto grid
template <typename detector_t, typename view_t>
auto material_grid(const detector_t& detector, const dindex index,
                   const view_t& view,
                   const styling::grid_style& style =
                       styling::tableau_colorblind::grid_style) {
  using scalar_t = dscalar<typename detector_t::algebra_type>;

  const auto& sf_desc = detector.surface(index);
  const auto& link = sf_desc.material();

  // Proactively calculate the reference radius for a cylinder grid
  // (will only be used if the volume actually holds a barrel grid)
  auto r =
      detector.mask_store().template visit<detray::detail::outer_radius_getter>(
          sf_desc.mask());

  scalar_t cyl_ref_radius{0.f};
  if (r.has_value()) {
    cyl_ref_radius = *r;
  }

  return svgtools::conversion::grid(detector.material_store(), link, view,
                                    cyl_ref_radius, style);
}

/// @brief A functor to fill the material proto grid with proto material slabs
struct material_converter {
  template <typename mat_coll_t, typename index_t,
            concepts::transform3D transform3_t>
  DETRAY_HOST inline auto operator()(const mat_coll_t& mat_coll,
                                     const index_t& index,
                                     const transform3_t&) const {
    using material_t = typename mat_coll_t::value_type;

    std::vector<std::vector<actsvg::proto::material_slab>> m_matrix;

    if constexpr (concepts::surface_material<material_t> &&
                  concepts::material_map<material_t>) {
      using loc_bin_idx_t = typename material_t::loc_bin_index;
      using algebra_t = typename material_t::local_frame_type::algebra_type;

      const auto material_map = mat_coll.at(index);

      // Create the bin associations
      dindex edges0 = material_map.template get_axis<0>().nbins();
      dindex edges1 = material_map.template get_axis<1>().nbins();

      // In the svg convention the phi axis has to be the second axis to
      // loop over
      constexpr bool is_cyl{
          std::is_same_v<typename material_t::local_frame_type,
                         detray::cylindrical2D<algebra_t>> ||
          std::is_same_v<typename material_t::local_frame_type,
                         detray::concentric_cylindrical2D<algebra_t>>};
      if constexpr (is_cyl) {
        dindex tmp = edges0;
        edges0 = edges1;
        edges1 = tmp;
      }

      m_matrix.reserve(edges1);

      // Material map is always 2-dimensional
      for (dindex j = 0u; j < edges1; ++j) {
        std::vector<actsvg::proto::material_slab> m_matrix_row;
        m_matrix_row.reserve(edges1);

        for (dindex i = 0u; i < edges0; ++i) {
          loc_bin_idx_t bin_idx{i, j};
          if constexpr (is_cyl) {
            bin_idx = {j, i};
          }

          const auto& mat_slab = material_map.bin(bin_idx).ref();
          const auto& mat = mat_slab.get_material();

          actsvg::proto::material_slab p_mat_slab;
          p_mat_slab._properties = {
              static_cast<actsvg::scalar>(mat.X0()),
              static_cast<actsvg::scalar>(mat.L0()),
              static_cast<actsvg::scalar>(mat.Ar()),
              static_cast<actsvg::scalar>(mat.Z()),
              static_cast<actsvg::scalar>(mat.mass_density()),
              static_cast<actsvg::scalar>(mat_slab.thickness())};

          m_matrix_row.push_back(std::move(p_mat_slab));
        }

        m_matrix.push_back(std::move(m_matrix_row));
      }
    }

    return m_matrix;
  }
};

/// @brief Calculate the proto surface material of a surface.
///
/// @param d_surface The detray surface.
/// @param context The geometry context.
///
/// @returns An actsvg proto surface material the material map.
template <typename detector_t, typename view_t>
auto surface_material(const detector_t& detector,
                      const detray::geometry::surface<detector_t>& d_surface,
                      const view_t& view,
                      const styling::surface_material_style& style =
                          styling::tableau_colorblind::material_style) {
  // Convert grid, if present
  auto [p_grid, grid_type] = svgtools::conversion::material_grid(
      detector, d_surface.index(), view, style._grid_style);

  actsvg::proto::surface_material p_material;

  // Create the surface material
  if (p_grid.has_value()) {
    // Fill material data
    auto m_matrix = d_surface.template visit_material<material_converter>(
        d_surface.transform({}));

    p_material = {m_matrix, *p_grid};
    p_material._material_ranges =
        actsvg::proto::material_ranges(p_material._material_matrix);

    svgtools::styling::apply_style(p_material, style);
  }

  return p_material;
}

}  // namespace detray::svgtools::conversion
