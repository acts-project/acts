// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/grid_axis.hpp"
#include "detray/geometry/coordinates/coordinates.hpp"
#include "detray/utils/grid/axis.hpp"
#include "detray/utils/grid/detail/axis_binning.hpp"
#include "detray/utils/grid/detail/axis_bounds.hpp"
#include "detray/utils/type_list.hpp"

// System include(s).
#include <type_traits>

namespace detray {

namespace axis::detail {

// Forward declarations of helper types
template <typename, axis::bounds, template <typename, typename> class,
          template <typename, typename> class,
          template <typename, typename> class>
struct get_axes_types;

template <typename, bool, typename, typename, typename = void>
struct get_coordinate_axes_type;

}  // namespace axis::detail

/// Template parameter for the @c grid: Get a the binning and bounds of grid
/// axes from geometric shapes
template <typename shape_t, axis::bounds e_bounds = axis::bounds::e_closed,
          template <typename, typename> class binning0_t = axis::regular,
          template <typename, typename> class binning1_t = axis::regular,
          template <typename, typename> class binning2_t = axis::regular>
struct axes {
  static constexpr auto dim{shape_t::dim};
  using bounds = void;
  template <concepts::algebra A>
  using type = axis::detail::get_axes_types<
      typename shape_t::template local_frame_type<A>, e_bounds, binning0_t,
      binning1_t, binning2_t>;
};

/// Helper trait to resolve the type of a @c multi_axis from a shape
template <typename axes_t, concepts::algebra algebra_t, bool is_owning = true,
          typename containers = host_container_types>
using coordinate_axes = typename axis::detail::get_coordinate_axes_type<
    axes_t, is_owning, containers, algebra_t>::type;

namespace axis::detail {

/// Determine axis bounds as either 'open' or 'closed' for non-circular axes.
template <axis::bounds s, axis::label axis_label>
using bounds_t =
    std::conditional_t<s == axis::bounds::e_open, axis::open<axis_label>,
                       axis::closed<axis_label>>;

/// Assemble a @c multi_axis type from a coordinate frame with given bounds and
/// binnings - unspecified type
template <typename frame_t, axis::bounds e_bounds,
          template <typename, typename> class binning0_t,
          template <typename, typename> class binning1_t,
          template <typename, typename> class binning2_t>
struct get_axes_types {};

/// Behaviour of the three local axes (linear in x and y)
template <concepts::algebra A, axis::bounds e_bounds,
          template <typename, typename> class binning0_t,
          template <typename, typename> class binning1_t,
          template <typename, typename> class binning2_t>
struct get_axes_types<cartesian2D<A>, e_bounds, binning0_t, binning1_t,
                      binning2_t> {
  template <concepts::algebra algebra_t>
  using frame = cartesian2D<algebra_t>;

  static constexpr auto label0{axis::label::e_x};
  static constexpr auto label1{axis::label::e_y};

  using bounds = detray::types::list<bounds_t<e_bounds, label0>,
                                     bounds_t<e_bounds, label1>>;
  template <concepts::scalar S, typename C>
  using binning = detray::types::list<binning0_t<S, C>, binning1_t<S, C>>;
};

/// Behaviour of the three local axes (linear in x, y, z)
template <concepts::algebra A, axis::bounds e_bounds,
          template <typename, typename> class binning0_t,
          template <typename, typename> class binning1_t,
          template <typename, typename> class binning2_t>
struct get_axes_types<cartesian3D<A>, e_bounds, binning0_t, binning1_t,
                      binning2_t> {
  template <concepts::algebra algebra_t>
  using frame = cartesian3D<algebra_t>;

  static constexpr auto label0{axis::label::e_x};
  static constexpr auto label1{axis::label::e_y};
  static constexpr auto label2{axis::label::e_z};

  using bounds = detray::types::list<bounds_t<e_bounds, label0>,
                                     bounds_t<e_bounds, label1>,
                                     bounds_t<e_bounds, label2>>;
  template <concepts::scalar S, typename C>
  using binning =
      detray::types::list<binning0_t<S, C>, binning1_t<S, C>, binning2_t<S, C>>;
};

/// Behaviour of the two local axes (linear in r, circular in phi)
template <concepts::algebra A, axis::bounds e_bounds,
          template <typename, typename> class binning0_t,
          template <typename, typename> class binning1_t,
          template <typename, typename> class binning2_t>
struct get_axes_types<polar2D<A>, e_bounds, binning0_t, binning1_t,
                      binning2_t> {
  template <concepts::algebra algebra_t>
  using frame = polar2D<algebra_t>;

  static constexpr auto label0{axis::label::e_r};
  static constexpr auto label1{axis::label::e_phi};

  using bounds =
      detray::types::list<bounds_t<e_bounds, label0>, axis::circular<label1>>;
  template <concepts::scalar S, typename C>
  using binning = detray::types::list<binning0_t<S, C>, binning1_t<S, C>>;
};

/// Behaviour of the two local axes (circular in r-phi, linear in z)
template <concepts::algebra A, axis::bounds e_bounds,
          template <typename, typename> class binning0_t,
          template <typename, typename> class binning1_t,
          template <typename, typename> class binning2_t>
struct get_axes_types<cylindrical2D<A>, e_bounds, binning0_t, binning1_t,
                      binning2_t> {
  template <concepts::algebra algebra_t>
  using frame = cylindrical2D<algebra_t>;

  static constexpr auto label0{axis::label::e_rphi};
  static constexpr auto label1{axis::label::e_cyl_z};

  using bounds =
      detray::types::list<axis::circular<label0>, bounds_t<e_bounds, label1>>;
  template <concepts::scalar S, typename C>
  using binning = detray::types::list<binning0_t<S, C>, binning1_t<S, C>>;
};

/// Behaviour of the two local axes (circular in r-phi, linear in z)
template <concepts::algebra A, axis::bounds e_bounds,
          template <typename, typename> class binning0_t,
          template <typename, typename> class binning1_t,
          template <typename, typename> class binning2_t>
struct get_axes_types<concentric_cylindrical2D<A>, e_bounds, binning0_t,
                      binning1_t, binning2_t> {
  template <concepts::algebra algebra_t>
  using frame = concentric_cylindrical2D<algebra_t>;

  static constexpr auto label0{axis::label::e_rphi};
  static constexpr auto label1{axis::label::e_cyl_z};

  using bounds =
      detray::types::list<axis::circular<label0>, bounds_t<e_bounds, label1>>;
  template <concepts::scalar S, typename C>
  using binning = detray::types::list<binning0_t<S, C>, binning1_t<S, C>>;
};

/// Behaviour of the two local axes (linear in r, circular in phi, linear in z)
template <concepts::algebra A, axis::bounds e_bounds,
          template <typename, typename> class binning0_t,
          template <typename, typename> class binning1_t,
          template <typename, typename> class binning2_t>
struct get_axes_types<cylindrical3D<A>, e_bounds, binning0_t, binning1_t,
                      binning2_t> {
  template <concepts::algebra algebra_t>
  using frame = cylindrical3D<algebra_t>;

  static constexpr auto label0{axis::label::e_r};
  static constexpr auto label1{axis::label::e_phi};
  static constexpr auto label2{axis::label::e_z};

  using bounds =
      detray::types::list<bounds_t<e_bounds, label0>, axis::circular<label1>,
                          bounds_t<e_bounds, label2>>;
  template <concepts::scalar S, typename C>
  using binning =
      detray::types::list<binning0_t<S, C>, binning1_t<S, C>, binning2_t<S, C>>;
};

/// Construct a @c multi_axis type from a given axes-shape
template <typename axes_t, bool is_owning, typename containers,
          concepts::algebra algebra_t>
  requires std::is_same_v<typename axes_t::bounds, void>
struct get_coordinate_axes_type<axes_t, is_owning, containers, algebra_t> {
  using type = typename axis::detail::multi_axis_assembler<
      is_owning, containers,
      typename axes_t::template type<algebra_t>::template frame<algebra_t>,
      typename axes_t::template type<algebra_t>::bounds,
      typename axes_t::template type<algebra_t>::template binning<
          dscalar<algebra_t>, containers>>::type;
};

/// Don't do anything if the type is already a @c multi_axis
template <typename axes_t, bool is_owning, typename containers,
          concepts::algebra algebra_t>
  requires std::is_object_v<typename axes_t::bounds>
struct get_coordinate_axes_type<axes_t, is_owning, containers, algebra_t> {
  using type = axes_t;
};

}  // namespace axis::detail

}  // namespace detray
