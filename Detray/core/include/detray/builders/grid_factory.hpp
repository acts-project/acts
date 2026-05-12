// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/utils/grid/axis.hpp"
#include "detray/utils/grid/concepts.hpp"
#include "detray/utils/grid/detail/axis_helpers.hpp"
#include "detray/utils/grid/grid.hpp"
#include "detray/utils/grid/grid_collection.hpp"
#include "detray/utils/grid/populators.hpp"
#include "detray/utils/grid/serializers.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <cassert>
#include <iostream>
#include <string>
#include <type_traits>

namespace detray {

/// @brief Provides functionality to instantiate grids and grid collections.
///
/// @tparam bin_t type bins in the grid
/// @tparam serialzier_t  type of the serializer to the storage representations
/// @tparam algebra_t the matrix/vector/point types to use
/// @tparam container_t the container types to use
///
/// throughout:
/// @note that bin_edges is only used for variable binning
/// @note that if non-zero axis_spans are provided the values of the
/// mask is overwritten
template <typename bin_t, template <std::size_t> class serializer_t,
          concepts::algebra algebra_t>
class grid_factory {
 public:
  // All grids are owning since they are used to fill the data
  static constexpr bool is_owning = true;

  using bin_type = bin_t;
  template <typename grid_shape_t>
  using grid_type = grid<algebra_t, axes<grid_shape_t>, bin_type, serializer_t>;
  template <typename grid_shape_t>
  using loc_bin_index = typename grid_type<grid_shape_t>::loc_bin_index;

  using scalar_type = dscalar<algebra_t>;
  template <typename T>
  using vector_type = host_container_types::template vector_type<T>;
  using algebra_type = algebra_t;

  grid_factory() = default;

  /// Takes the resource of the detector to allocate memory correctly
  explicit grid_factory(vecmem::memory_resource &resource)
      : m_resource(&resource) {}

  //
  // annulus 2D
  //
  template <
      typename r_bounds = axis::closed<axis::label::e_r>,
      typename phi_bounds = axis::circular<>,
      typename r_binning = axis::regular<scalar_type, host_container_types>,
      typename phi_binning = axis::regular<scalar_type, host_container_types>>
    requires std::is_enum_v<decltype(r_bounds::label)>
  auto new_grid(const mask<annulus2D, algebra_type> &grid_bounds,
                const darray<std::size_t, 2UL> n_bins,
                const std::vector<std::pair<loc_bin_index<annulus2D>, dindex>>
                    &bin_capacities = {},
                const darray<std::vector<scalar_type>, 2UL> &bin_edges =
                    darray<std::vector<scalar_type>, 2UL>(),
                const darray<std::vector<scalar_type>, 2UL> &axis_spans =
                    darray<std::vector<scalar_type>, 2UL>()) const {
    static_assert(
        std::is_same_v<phi_bounds, axis::circular<>>,
        "Phi axis bounds need to be circular for stereo annulus shape");

    // Axes boundaries and local indices
    using boundary = annulus2D::boundaries;
    using axes_t = axes<annulus2D>::template type<algebra_t>;
    using local_frame = typename axes_t::template frame<algebra_t>;

    constexpr auto e_r_axis = static_cast<dindex>(axes_t::label0);
    constexpr auto e_phi_axis = static_cast<dindex>(axes_t::label1);

    auto b_values = grid_bounds.values();
    // Overwrite the mask values if axis spans are provided
    if (!axis_spans[0UL].empty()) {
      assert(axis_spans[0UL].size() == 2UL);
      b_values[boundary::e_min_r] = axis_spans[0UL].at(0UL);
      b_values[boundary::e_max_r] = axis_spans[0UL].at(1UL);
    }
    scalar_type min_phi =
        b_values[boundary::e_average_phi] - b_values[boundary::e_min_phi_rel];
    scalar_type max_phi =
        b_values[boundary::e_average_phi] + b_values[boundary::e_max_phi_rel];
    if (!axis_spans[1UL].empty()) {
      assert(axis_spans[1UL].size() == 2UL);
      min_phi = axis_spans[1UL].at(0UL);
      max_phi = axis_spans[1UL].at(1UL);
    }

    return new_grid<local_frame>(
        {b_values[boundary::e_min_r], b_values[boundary::e_max_r], min_phi,
         max_phi},
        {n_bins[e_r_axis], n_bins[e_phi_axis]}, bin_capacities,
        {bin_edges[e_r_axis], bin_edges[e_phi_axis]},
        types::list<r_bounds, phi_bounds>{},
        types::list<r_binning, phi_binning>{});
  }

  //
  // cuboid 3D
  //
  template <
      typename x_bounds = axis::closed<axis::label::e_x>,
      typename y_bounds = axis::closed<axis::label::e_y>,
      typename z_bounds = axis::closed<axis::label::e_z>,
      typename x_binning = axis::regular<scalar_type, host_container_types>,
      typename y_binning = axis::regular<scalar_type, host_container_types>,
      typename z_binning = axis::regular<scalar_type, host_container_types>>
    requires std::is_enum_v<decltype(x_bounds::label)>
  auto new_grid(const mask<cuboid3D, algebra_type> &grid_bounds,
                const darray<std::size_t, 3UL> n_bins,
                const std::vector<std::pair<loc_bin_index<cuboid3D>, dindex>>
                    &bin_capacities = {},
                const darray<std::vector<scalar_type>, 3UL> &bin_edges =
                    darray<std::vector<scalar_type>, 3UL>(),
                const darray<std::vector<scalar_type>, 3UL> &axis_spans =
                    darray<std::vector<scalar_type>, 3UL>()) const {
    // Axes boundaries and local indices
    using boundary = cuboid3D::boundaries;
    using axes_t = axes<cuboid3D>::template type<algebra_t>;
    using local_frame = typename axes_t::template frame<algebra_t>;

    constexpr auto e_x_axis = static_cast<dindex>(axes_t::label0);
    constexpr auto e_y_axis = static_cast<dindex>(axes_t::label1);
    constexpr auto e_z_axis = static_cast<dindex>(axes_t::label2);

    // Overwrite the mask values if axis spans are provided
    auto b_values = grid_bounds.values();
    if (!axis_spans[0UL].empty()) {
      assert(axis_spans[0UL].size() == 2UL);
      b_values[boundary::e_min_x] = axis_spans[0UL].at(0UL);
      b_values[boundary::e_max_x] = axis_spans[0UL].at(1UL);
    }
    if (!axis_spans[1UL].empty()) {
      assert(axis_spans[1UL].size() == 2UL);
      b_values[boundary::e_min_y] = axis_spans[1UL].at(0UL);
      b_values[boundary::e_max_y] = axis_spans[1UL].at(1UL);
    }
    if (!axis_spans[2UL].empty()) {
      assert(axis_spans[2UL].size() == 2UL);
      b_values[boundary::e_min_z] = axis_spans[2UL].at(0UL);
      b_values[boundary::e_max_z] = axis_spans[2UL].at(1UL);
    }

    return new_grid<local_frame>(
        {b_values[boundary::e_min_x], b_values[boundary::e_max_x],
         b_values[boundary::e_min_y], b_values[boundary::e_max_y],
         b_values[boundary::e_min_z], b_values[boundary::e_max_z]},
        {n_bins[e_x_axis], n_bins[e_y_axis], n_bins[e_z_axis]}, bin_capacities,
        {bin_edges[e_x_axis], bin_edges[e_y_axis], bin_edges[e_z_axis]},
        types::list<x_bounds, y_bounds, z_bounds>{},
        types::list<x_binning, y_binning, z_binning>{});
  }

  //
  // cylinder 2D
  //
  template <
      typename rphi_bounds = axis::circular<axis::label::e_rphi>,
      typename z_bounds = axis::closed<axis::label::e_cyl_z>,
      typename rphi_binning = axis::regular<scalar_type, host_container_types>,
      typename z_binning = axis::regular<scalar_type, host_container_types>>
    requires std::is_enum_v<decltype(rphi_bounds::label)>
  auto new_grid(const mask<cylinder2D, algebra_type> &grid_bounds,
                const darray<std::size_t, 2UL> n_bins,
                const std::vector<std::pair<loc_bin_index<cylinder2D>, dindex>>
                    &bin_capacities = {},
                const darray<std::vector<scalar_type>, 2UL> &bin_edges =
                    darray<std::vector<scalar_type>, 2UL>(),
                const darray<std::vector<scalar_type>, 2UL> &axis_spans =
                    darray<std::vector<scalar_type>, 2UL>()) const {
    static_assert(
        std::is_same_v<rphi_bounds, axis::circular<axis::label::e_rphi>>,
        "Phi axis bounds need to be circular for cylinder2D shape");

    // Axes boundaries and local indices
    using boundary = cylinder2D::boundaries;
    using axes_t = axes<cylinder2D>::template type<algebra_t>;
    using local_frame = typename axes_t::template frame<algebra_t>;

    constexpr auto e_rphi_axis = static_cast<dindex>(axes_t::label0);
    constexpr auto e_z_axis = static_cast<dindex>(axes_t::label1);

    auto b_values = grid_bounds.values();
    if (!axis_spans[1UL].empty()) {
      assert(axis_spans[1UL].size() == 2UL);
      b_values[boundary::e_lower_z] = axis_spans[1UL].at(0UL);
      b_values[boundary::e_upper_z] = axis_spans[1UL].at(1UL);
    }
    return new_grid<local_frame>(
        {-constant<scalar_type>::pi, constant<scalar_type>::pi,
         b_values[boundary::e_lower_z], b_values[boundary::e_upper_z]},
        {n_bins[e_rphi_axis], n_bins[e_z_axis]}, bin_capacities,
        {bin_edges[e_rphi_axis], bin_edges[e_z_axis]},
        types::list<rphi_bounds, z_bounds>{},
        types::list<rphi_binning, z_binning>{});
  }

  //
  // concentric cylinder 2D
  //
  template <
      typename rphi_bounds = axis::circular<axis::label::e_rphi>,
      typename z_bounds = axis::closed<axis::label::e_cyl_z>,
      typename rphi_binning = axis::regular<scalar_type, host_container_types>,
      typename z_binning = axis::regular<scalar_type, host_container_types>>
    requires std::is_enum_v<decltype(rphi_bounds::label)>
  auto new_grid(
      const mask<concentric_cylinder2D, algebra_type> &grid_bounds,
      const darray<std::size_t, 2UL> n_bins,
      const std::vector<std::pair<loc_bin_index<concentric_cylinder2D>, dindex>>
          &bin_capacities = {},
      const darray<std::vector<scalar_type>, 2UL> &bin_edges =
          darray<std::vector<scalar_type>, 2UL>(),
      const darray<std::vector<scalar_type>, 2UL> &axis_spans =
          darray<std::vector<scalar_type>, 2UL>()) const {
    static_assert(
        std::is_same_v<rphi_bounds, axis::circular<axis::label::e_rphi>>,
        "Phi axis bounds need to be circular for cylinder2D portal shape");

    // Axes boundaries and local indices
    using boundary = concentric_cylinder2D::boundaries;
    using axes_t = axes<concentric_cylinder2D>::template type<algebra_t>;
    using local_frame = typename axes_t::template frame<algebra_t>;

    constexpr auto e_rphi_axis = static_cast<dindex>(axes_t::label0);
    constexpr auto e_z_axis = static_cast<dindex>(axes_t::label1);

    auto b_values = grid_bounds.values();
    if (!axis_spans[1UL].empty()) {
      assert(axis_spans[1UL].size() == 2UL);
      b_values[boundary::e_lower_z] = axis_spans[1UL].at(0UL);
      b_values[boundary::e_upper_z] = axis_spans[1UL].at(1UL);
    }

    return new_grid<local_frame>(
        {-constant<scalar_type>::pi, constant<scalar_type>::pi,
         b_values[boundary::e_lower_z], b_values[boundary::e_upper_z]},
        {n_bins[e_rphi_axis], n_bins[e_z_axis]}, bin_capacities,
        {bin_edges[e_rphi_axis], bin_edges[e_z_axis]},
        types::list<rphi_bounds, z_bounds>{},
        types::list<rphi_binning, z_binning>{});
  }

  //
  // cylinder 3D
  //
  template <
      typename r_bounds = axis::closed<axis::label::e_r>,
      typename phi_bounds = axis::circular<>,
      typename z_bounds = axis::closed<axis::label::e_z>,
      typename r_binning = axis::regular<scalar_type, host_container_types>,
      typename phi_binning = axis::regular<scalar_type, host_container_types>,
      typename z_binning = axis::regular<scalar_type, host_container_types>>
    requires std::is_enum_v<decltype(r_bounds::label)>
  auto new_grid(const mask<cylinder3D, algebra_type> &grid_bounds,
                const darray<std::size_t, 3UL> n_bins,
                const std::vector<std::pair<loc_bin_index<cylinder3D>, dindex>>
                    &bin_capacities = {},
                const darray<std::vector<scalar_type>, 3UL> &bin_edges =
                    darray<std::vector<scalar_type>, 3UL>(),
                const darray<std::vector<scalar_type>, 3UL> &axis_spans =
                    darray<std::vector<scalar_type>, 3UL>()) const {
    static_assert(std::is_same_v<phi_bounds, axis::circular<>>,
                  "Phi axis bounds need to be circular for cylinder3D shape");

    // Axes boundaries and local indices
    using boundary = cylinder3D::boundaries;
    using axes_t = axes<cylinder3D>::template type<algebra_t>;
    using local_frame = typename axes_t::template frame<algebra_t>;

    constexpr auto e_r_axis = static_cast<dindex>(axes_t::label0);
    constexpr auto e_phi_axis = static_cast<dindex>(axes_t::label1);
    constexpr auto e_z_axis = static_cast<dindex>(axes_t::label2);

    auto b_values = grid_bounds.values();

    // Overwrite the mask values if axis spans are provided
    if (!axis_spans[0UL].empty()) {
      assert(axis_spans[0UL].size() == 2UL);
      b_values[boundary::e_min_r] = axis_spans[0UL].at(0UL);
      b_values[boundary::e_max_r] = axis_spans[0UL].at(1UL);
    }
    scalar_type min_phi = -constant<scalar_type>::pi;
    scalar_type max_phi = constant<scalar_type>::pi;
    if (!axis_spans[1UL].empty()) {
      assert(axis_spans[1UL].size() == 2UL);
      min_phi = axis_spans[1UL].at(0UL);
      max_phi = axis_spans[1UL].at(1UL);
    }
    if (!axis_spans[2UL].empty()) {
      assert(axis_spans[2UL].size() == 2UL);
      b_values[boundary::e_min_z] = axis_spans[2UL].at(0UL);
      b_values[boundary::e_max_z] = axis_spans[2UL].at(1UL);
    }

    return new_grid<local_frame>(
        {b_values[boundary::e_min_r], b_values[boundary::e_max_r], min_phi,
         max_phi, -b_values[boundary::e_min_z], b_values[boundary::e_max_z]},
        {n_bins[e_r_axis], n_bins[e_phi_axis], n_bins[e_z_axis]},
        bin_capacities,
        {bin_edges[e_r_axis], bin_edges[e_phi_axis], bin_edges[e_z_axis]},
        types::list<r_bounds, phi_bounds, z_bounds>{},
        types::list<r_binning, phi_binning, z_binning>{});
  }

  //
  // polar 2D
  //
  template <
      typename r_bounds = axis::closed<axis::label::e_r>,
      typename phi_bounds = axis::circular<>,
      typename r_binning = axis::regular<scalar_type, host_container_types>,
      typename phi_binning = axis::regular<scalar_type, host_container_types>>
    requires std::is_enum_v<decltype(r_bounds::label)>
  auto new_grid(const mask<ring2D, algebra_type> &grid_bounds,
                const darray<std::size_t, 2UL> n_bins,
                const std::vector<std::pair<loc_bin_index<ring2D>, dindex>>
                    &bin_capacities = {},
                const darray<std::vector<scalar_type>, 2UL> &bin_edges =
                    darray<std::vector<scalar_type>, 2UL>(),
                const darray<std::vector<scalar_type>, 2UL> &axis_spans =
                    darray<std::vector<scalar_type>, 2UL>()) const {
    static_assert(std::is_same_v<phi_bounds, axis::circular<>>,
                  "Phi axis bounds need to be circular for ring shape");

    // Axes boundaries and local indices
    using boundary = ring2D::boundaries;
    using axes_t = axes<ring2D>::template type<algebra_t>;
    using local_frame = typename axes_t::template frame<algebra_t>;

    constexpr auto e_r_axis = static_cast<dindex>(axes_t::label0);
    constexpr auto e_phi_axis = static_cast<dindex>(axes_t::label1);

    auto b_values = grid_bounds.values();
    // Overwrite the mask values if axis spans are provided
    if (!axis_spans[0UL].empty()) {
      assert(axis_spans[0UL].size() == 2UL);
      b_values[boundary::e_inner_r] = axis_spans[0UL].at(0UL);
      b_values[boundary::e_outer_r] = axis_spans[0UL].at(1UL);
    }

    return new_grid<local_frame>(
        {b_values[boundary::e_inner_r], b_values[boundary::e_outer_r],
         -constant<scalar_type>::pi, constant<scalar_type>::pi},
        {n_bins[e_r_axis], n_bins[e_phi_axis]}, bin_capacities,
        {bin_edges[e_r_axis], bin_edges[e_phi_axis]},
        types::list<r_bounds, phi_bounds>{},
        types::list<r_binning, phi_binning>{});
  }

  //
  // rectangle 2D
  //
  template <
      typename x_bounds = axis::closed<axis::label::e_x>,
      typename y_bounds = axis::closed<axis::label::e_y>,
      typename x_binning = axis::regular<scalar_type, host_container_types>,
      typename y_binning = axis::regular<scalar_type, host_container_types>>
    requires std::is_enum_v<decltype(x_bounds::label)>
  auto new_grid(const mask<rectangle2D, algebra_type> &grid_bounds,
                const darray<std::size_t, 2UL> n_bins,
                const std::vector<std::pair<loc_bin_index<rectangle2D>, dindex>>
                    &bin_capacities = {},
                const darray<std::vector<scalar_type>, 2UL> &bin_edges =
                    darray<std::vector<scalar_type>, 2UL>(),
                const darray<std::vector<scalar_type>, 2UL> &axis_spans =
                    darray<std::vector<scalar_type>, 2UL>()) const {
    // Axes boundaries and local indices
    using boundary = rectangle2D::boundaries;
    using axes_t = axes<rectangle2D>::template type<algebra_t>;
    using local_frame = typename axes_t::template frame<algebra_t>;

    constexpr auto e_x_axis = static_cast<dindex>(axes_t::label0);
    constexpr auto e_y_axis = static_cast<dindex>(axes_t::label1);

    auto b_values = grid_bounds.values();
    // Overwrite the mask values if axis spans are provided
    if (!axis_spans[0UL].empty()) {
      assert(axis_spans[0UL].size() == 2UL);
      b_values[boundary::e_half_x] = axis_spans[0UL].at(1UL);
    }
    if (!axis_spans[1UL].empty()) {
      assert(axis_spans[1UL].size() == 2UL);
      b_values[boundary::e_half_y] = axis_spans[1UL].at(1UL);
    }

    return new_grid<local_frame>(
        {-b_values[boundary::e_half_x], b_values[boundary::e_half_x],
         -b_values[boundary::e_half_y], b_values[boundary::e_half_y]},
        {n_bins[e_x_axis], n_bins[e_y_axis]}, bin_capacities,
        {bin_edges[e_x_axis], bin_edges[e_y_axis]},
        types::list<x_bounds, y_bounds>{}, types::list<x_binning, y_binning>{});
  }

  //
  // trapezoid 2D
  //
  template <
      typename x_bounds = axis::closed<axis::label::e_x>,
      typename y_bounds = axis::closed<axis::label::e_y>,
      typename x_binning = axis::regular<scalar_type, host_container_types>,
      typename y_binning = axis::regular<scalar_type, host_container_types>>
    requires std::is_enum_v<decltype(x_bounds::label)>
  auto new_grid(const mask<trapezoid2D, algebra_type> &grid_bounds,
                const darray<std::size_t, 2UL> n_bins,
                const std::vector<std::pair<loc_bin_index<trapezoid2D>, dindex>>
                    &bin_capacities = {},
                const darray<std::vector<scalar_type>, 2UL> &bin_edges =
                    darray<std::vector<scalar_type>, 2UL>(),
                const darray<std::vector<scalar_type>, 2UL> &axis_spans =
                    darray<std::vector<scalar_type>, 2UL>()) const {
    // Axes boundaries and local indices
    using boundary = trapezoid2D::boundaries;
    using axes_t = axes<trapezoid2D>::template type<algebra_t>;
    using local_frame = typename axes_t::template frame<algebra_t>;

    constexpr auto e_x_axis = static_cast<dindex>(axes_t::label0);
    constexpr auto e_y_axis = static_cast<dindex>(axes_t::label1);

    auto b_values = grid_bounds.values();
    // Overwrite the mask values if axis spans are provided
    if (!axis_spans[0UL].empty()) {
      assert(axis_spans[0UL].size() == 2UL);
      b_values[boundary::e_half_length_1] = axis_spans[0UL].at(1UL);
    }
    if (!axis_spans[1UL].empty()) {
      assert(axis_spans[1UL].size() == 2UL);
      b_values[boundary::e_half_length_2] = axis_spans[1UL].at(1UL);
    }

    return new_grid<local_frame>(
        {-b_values[boundary::e_half_length_1],
         b_values[boundary::e_half_length_1],
         -b_values[boundary::e_half_length_2],
         b_values[boundary::e_half_length_2]},
        {n_bins[e_x_axis], n_bins[e_y_axis]}, bin_capacities,
        {bin_edges[e_x_axis], bin_edges[e_y_axis]},
        types::list<x_bounds, y_bounds>{}, types::list<x_binning, y_binning>{});
  }

  /// @brief Create and empty grid with fully initialized axes.
  ///
  /// @tparam grid_shape_t the shape of the resulting grid
  ///         (e.g. cylinder2D).
  /// @tparam e_bounds the bounds of the regular axes
  ///         (open vs. closed bounds).
  /// @tparam binnings the binning types of the axes
  ///         (regular vs. irregular)
  ///
  /// @param spans the span of the axis values for regular axes, otherwise
  ///              ignored.
  /// @param n_bins the number of bins for regular axes, otherwise ignored
  /// @param ax_bin_edges the explicit bin edges for irregular axes
  ///                     (lower bin edges + the the upper edge of the
  ///                     last bin), otherwise ignored.
  template <concepts::coordinate_frame grid_frame_t, typename... bound_ts,
            typename... binning_ts>
  auto new_grid(
      const std::vector<scalar_type> &spans,
      const std::vector<std::size_t> &n_bins,
      const std::vector<std::pair<axis::multi_bin<sizeof...(bound_ts)>, dindex>>
          &bin_capacities = {},
      const std::vector<std::vector<scalar_type>> &ax_bin_edges = {},
      const types::list<bound_ts...> & /*unused*/ = {},
      const types::list<binning_ts...> & /*unused*/ = {}) const {
    static_assert(sizeof...(bound_ts) == sizeof...(binning_ts),
                  "number of axis bounds and binning types has to match");

    // Build the coordinate axes and the grid
    using axes_t = axis::multi_axis<is_owning, grid_frame_t,
                                    axis::single_axis<bound_ts, binning_ts>...>;
    using grid_t = grid<algebra_t, axes_t, bin_t, serializer_t>;

    return new_grid<grid_t>(spans, n_bins, bin_capacities, ax_bin_edges);
  }

  /// Helper to build grid from shape plus binning and bounds types
  template <typename grid_shape_t, typename... bound_ts, typename... binning_ts>
    requires concepts::shape<grid_shape_t, algebra_t>
  auto new_grid(
      const std::vector<scalar_type> &spans,
      const std::vector<std::size_t> &n_bins,
      const std::vector<std::pair<axis::multi_bin<sizeof...(bound_ts)>, dindex>>
          &bin_capacities = {},
      const std::vector<std::vector<scalar_type>> &ax_bin_edges = {},
      const types::list<bound_ts...> &bounds = {},
      const types::list<binning_ts...> &binnings = {}) const {
    return new_grid<
        typename grid_shape_t::template local_frame_type<algebra_t>>(
        spans, n_bins, bin_capacities, ax_bin_edges, bounds, binnings);
  }

  /// Helper overload for grid builder: Build from mask and resolve bounds
  /// and binnings from concrete grid type
  template <concepts::grid grid_t, typename grid_shape_t>
    requires concepts::shape<grid_shape_t, algebra_t>
  auto new_grid(
      const mask<grid_shape_t, algebra_type> &m,
      const darray<std::size_t, grid_t::dim> &n_bins,
      const std::vector<std::pair<typename grid_t::loc_bin_index, dindex>>
          &bin_capacities = {},
      const darray<std::vector<scalar_type>, grid_t::dim> &ax_bin_edges = {})
      const {
    return new_grid(m, n_bins, bin_capacities, ax_bin_edges,
                    typename grid_t::axes_type::bounds{},
                    typename grid_t::axes_type::binnings{});
  }

  /// Helper overload for grid builder: Build from mask and resolve bounds
  /// and binnings
  template <typename grid_shape_t, typename... bound_ts, typename... binning_ts>
    requires concepts::shape<grid_shape_t, algebra_t>
  auto new_grid(
      const mask<grid_shape_t, algebra_type> &m,
      const darray<std::size_t, grid_shape_t::dim> &n_bins,
      const std::vector<std::pair<axis::multi_bin<sizeof...(bound_ts)>, dindex>>
          &bin_capacities = {},
      const darray<std::vector<scalar_type>, grid_shape_t::dim> &ax_bin_edges =
          {},
      const types::list<bound_ts...> & /*unused*/ = {},
      const types::list<binning_ts...> & /*unused*/ = {}) const {
    return new_grid<bound_ts..., binning_ts...>(m, n_bins, bin_capacities,
                                                ax_bin_edges);
  }

  /// @brief Create and empty grid with fully initialized axes.
  ///
  /// @tparam grid_t the type of the resulting grid
  ///
  /// @param spans the span of the axis values for regular axes, otherwise
  ///              ignored.
  /// @param n_bins the number of bins for regular axes, otherwise ignored
  /// @param ax_bin_edges the explicit bin edges for irregular axes
  ///                     (lower bin edges + the the upper edge of the
  ///                     last bin), otherwise ignored.
  template <concepts::grid grid_t>
  auto new_grid(
      const std::vector<scalar_type> &spans,
      const std::vector<std::size_t> &n_bins,
      [[maybe_unused]] const std::vector<
          std::pair<typename grid_t::loc_bin_index, dindex>> &bin_capacities =
          {},
      const std::vector<std::vector<scalar_type>> &ax_bin_edges = {}) const {
    using owning_grid_t = typename grid_t::template type<true>;
    using axes_t = typename owning_grid_t::axes_type;

    // Prepare data
    vector_type<dsized_index_range> axes_data{};
    vector_type<scalar_type> bin_edges{};

    // Call init for every axis
    unroll_axis_init<typename axes_t::binnings>(
        spans, n_bins, ax_bin_edges, axes_data, bin_edges,
        std::make_index_sequence<axes_t::dim>{});

    // Assemble the grid and return it
    axes_t axes(std::move(axes_data), std::move(bin_edges));

    typename owning_grid_t::bin_container_type bin_data{};

    // Bins with dynamic capacity and "index grids" need different treatment
    if constexpr (std::is_same_v<bin_t, bins::dynamic_array<
                                            typename grid_t::value_type>>) {
      // Bin and bin entries vector are separate containers
      bin_data.bins.resize(axes.nbins());
      // Set correct bin capacities, if they were provided
      if (!bin_capacities.empty()) {
        // Set memory offsets correctly
        typename grid_t::template serializer_type<grid_t::dim> serializer{};
        dindex total_cap{0u};
        for (const auto &capacity : bin_capacities) {
          // Get the empty bin by its global index
          auto &data = bin_data.bins[serializer(axes, capacity.first)];

          data.offset = total_cap;
          data.capacity = capacity.second;
          data.size = 0u;

          total_cap += data.capacity;
        }
        assert(total_cap > 0);
        bin_data.entries.resize(
            total_cap, detail::invalid_value<typename grid_t::value_type>());
      }
    } else {
      // Bin entries are contained in the bins directly
      bin_data.resize(axes.nbins(), bin_type{});
    }

    return owning_grid_t(std::move(bin_data), std::move(axes));
  }

 private:
  /// Initialize a single axis (either regular or irregular)
  /// @note change to template lambda as soon as it becomes available.
  template <std::size_t I, typename binnings>
  auto axis_init([[maybe_unused]] const std::vector<scalar_type> &spans,
                 [[maybe_unused]] const std::vector<std::size_t> &n_bins,
                 [[maybe_unused]] const std::vector<std::vector<scalar_type>>
                     &ax_bin_edges,
                 vector_type<dsized_index_range> &axes_data,
                 vector_type<scalar_type> &bin_edges) const {
    if constexpr (std::is_same_v<
                      types::at<binnings, I>,
                      axis::regular<scalar_type, host_container_types>>) {
      axes_data.push_back({static_cast<dindex>(bin_edges.size()),
                           static_cast<dindex>(n_bins.at(I))});
      bin_edges.push_back(spans.at(I * 2u));
      bin_edges.push_back(spans.at(I * 2u + 1u));
    } else {
      const auto &bin_edges_loc = ax_bin_edges.at(I);
      axes_data.push_back({static_cast<dindex>(bin_edges.size()),
                           static_cast<dindex>(bin_edges_loc.size() - 1u)});
      bin_edges.insert(bin_edges.end(), bin_edges_loc.begin(),
                       bin_edges_loc.end());
    }
  }

  /// Call axis init for every dimension
  template <typename binnings, std::size_t... I>
  auto unroll_axis_init(
      const std::vector<scalar_type> &spans,
      const std::vector<std::size_t> &n_bins,
      const std::vector<std::vector<scalar_type>> &ax_bin_edges,
      vector_type<dsized_index_range> &axes_data,
      vector_type<scalar_type> &bin_edges,
      std::index_sequence<I...> /*ids*/) const {
    (axis_init<I, binnings>(spans, n_bins, ax_bin_edges, axes_data, bin_edges),
     ...);
  }

  vecmem::memory_resource *m_resource{};
};

// Infer a grid factory type from an already completely assembled grid type
template <concepts::grid grid_t>
using grid_factory_type =
    grid_factory<typename grid_t::bin_type,
                 simple_serializer /*grid_t::template serializer_type*/,
                 typename grid_t::algebra_type>;

}  // namespace detray
