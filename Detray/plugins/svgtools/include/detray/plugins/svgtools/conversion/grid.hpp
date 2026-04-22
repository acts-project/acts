// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/grid_axis.hpp"
#include "detray/definitions/units.hpp"
#include "detray/utils/grid/concepts.hpp"

// Plugin include(s)
#include "detray/plugins/svgtools/styling/styling.hpp"
#include "detray/plugins/svgtools/utils/surface_kernels.hpp"

// Actsvg include(s)
#include "actsvg/proto/grid.hpp"

// System include(s)
#include <algorithm>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace detray::svgtools::conversion {

namespace detail {

enum class grid_type : std::uint_least8_t {
  e_barrel = 0,
  e_endcap = 1,
  e_unknown = 2
};

/// @returns the actsvg grid type and edge values for a detray 2D cylinder grid.
template <concepts::grid grid_t, typename view_t>
  requires std::is_same_v<
               typename grid_t::local_frame_type,
               detray::concentric_cylindrical2D<
                   typename grid_t::local_frame_type::algebra_type>> ||
           std::is_same_v<typename grid_t::local_frame_type,
                          detray::cylindrical2D<
                              typename grid_t::local_frame_type::algebra_type>>
inline auto grid_type_and_edges(const grid_t& grid, const view_t&) {
  using scalar_t = typename grid_t::local_frame_type::scalar_type;
  using axis_label = detray::axis::label;

  auto edges_phi = grid.template get_axis<axis_label::e_rphi>().bin_edges();
  auto edges_z = grid.template get_axis<axis_label::e_cyl_z>().bin_edges();
  // Unknown for 2D cylinder
  dvector<scalar_t> edges_r{};

  if constexpr (std::is_same_v<view_t, actsvg::views::x_y>) {
    return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_r_phi,
                      edges_r, edges_phi);
  }
  if constexpr (std::is_same_v<view_t, actsvg::views::z_r>) {
    return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_x_y, edges_z,
                      edges_r);
  }
  if constexpr (std::is_same_v<view_t, actsvg::views::z_phi> ||
                std::is_same_v<view_t, typename actsvg::views::z_rphi>) {
    return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_z_phi,
                      edges_z, edges_phi);
  }

  return std::tuple(grid_type::e_barrel, actsvg::proto::grid::e_x_y,
                    dvector<scalar_t>{}, dvector<scalar_t>{});
}

/// @returns the actsvg grid type and edge values for a detray disc grid.
template <concepts::grid grid_t, typename view_t>
  requires std::is_same_v<
      typename grid_t::local_frame_type,
      detray::polar2D<typename grid_t::local_frame_type::algebra_type>>
inline auto grid_type_and_edges(const grid_t& grid, const view_t&) {
  using scalar_t = typename grid_t::local_frame_type::scalar_type;
  using axis_label = detray::axis::label;

  auto edges_r = grid.template get_axis<axis_label::e_r>().bin_edges();
  auto edges_phi = grid.template get_axis<axis_label::e_phi>().bin_edges();

  if constexpr (std::is_same_v<view_t, typename actsvg::views::x_y>) {
    return std::tuple(grid_type::e_endcap, actsvg::proto::grid::e_r_phi,
                      edges_r, edges_phi);
  }

  return std::tuple(grid_type::e_endcap, actsvg::proto::grid::e_x_y, edges_r,
                    dvector<scalar_t>{});
}

/// @returns the actsvg grid type and edge values for a detray rectangular grid.
template <concepts::grid grid_t, typename view_t>
  requires std::is_same_v<
      typename grid_t::local_frame_type,
      detray::cartesian2D<typename grid_t::local_frame_type::algebra_type>>
inline auto grid_type_and_edges(const grid_t& grid, const view_t&) {
  using scalar_t = typename grid_t::local_frame_type::scalar_type;
  using axis_label = detray::axis::label;

  auto edges_x = grid.template get_axis<axis_label::e_x>().bin_edges();
  auto edges_y = grid.template get_axis<axis_label::e_y>().bin_edges();

  if constexpr (std::is_same_v<view_t, typename actsvg::views::x_y>) {
    return std::tuple(grid_type::e_endcap, actsvg::proto::grid::e_x_y, edges_x,
                      edges_y);
  }

  return std::tuple(grid_type::e_endcap, actsvg::proto::grid::e_x_y, edges_y,
                    dvector<scalar_t>{});
}

/// A functor to access the type and bin edges of a grid.
template <concepts::scalar scalar_t>
struct type_and_edge_getter {
  template <typename group_t, typename index_t, typename view_t>
  DETRAY_HOST_DEVICE inline auto operator()(
      [[maybe_unused]] const group_t& group,
      [[maybe_unused]] const index_t index,
      [[maybe_unused]] const view_t& view) const {
    using value_t = typename group_t::value_type;

    if constexpr (concepts::grid<value_t>) {
      // Only two dimensional grids for actsvg
      if constexpr (value_t::dim == 2) {
        return grid_type_and_edges(group.at(index), view);
      }
    }

    return std::tuple(grid_type::e_unknown, actsvg::proto::grid::e_x_y,
                      dvector<scalar_t>{}, dvector<scalar_t>{});
  }
};

}  // namespace detail

/// @brief Converts a detray grid to a actsvg proto grid.
///
/// @param store the data store that contains the grid
/// @param link the type id and index of the grid to be converted
/// @param view the view
/// @param ref_radius radius of the grid (only needed for cylindrical grids)
/// @param style the style settings
///
/// @returns a proto grid
template <typename store_t, typename link_t, typename view_t,
          concepts::scalar scalar_t>
auto grid(const store_t& store, const link_t& link, const view_t& view,
          const scalar_t ref_radius,
          const styling::grid_style& style =
              styling::tableau_colorblind::grid_style) {
  if (link.is_invalid()) {
    return std::tuple(std::optional<actsvg::proto::grid>{},
                      detail::grid_type::e_unknown);
  }

  auto [gr_type, view_type, edges0, edges1] =
      store.template visit<detail::type_and_edge_getter<scalar_t>>(link, view);
  actsvg::proto::grid p_grid;
  p_grid._type = view_type;

  // Find the correct grid radius
  if (gr_type == detail::grid_type::e_barrel) {
    p_grid._reference_r = static_cast<actsvg::scalar>(ref_radius);

    // Add the cylinder radius to the axis binning
    if constexpr (std::is_same_v<view_t, actsvg::views::x_y>) {
      if (edges0.empty()) {
        edges0 = {p_grid._reference_r, p_grid._reference_r};
      }
    }
    if constexpr (std::is_same_v<view_t, actsvg::views::z_r>) {
      if (edges1.empty()) {
        edges1 = {p_grid._reference_r, p_grid._reference_r};
      }
    }

  } else if (gr_type == detail::grid_type::e_endcap) {
    // An axis is always sorted
    p_grid._reference_r = static_cast<actsvg::scalar>(edges0.back());
  }

  std::ranges::transform(
      edges0, std::back_inserter(p_grid._edges_0),
      [](scalar_t v) { return static_cast<actsvg::scalar>(v); });
  std::ranges::transform(
      edges1, std::back_inserter(p_grid._edges_1),
      [](scalar_t v) { return static_cast<actsvg::scalar>(v); });

  svgtools::styling::apply_style(p_grid, style);

  return std::tuple(std::optional<actsvg::proto::grid>{p_grid}, gr_type);
}

}  // namespace detray::svgtools::conversion
