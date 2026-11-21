// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/array_soa_approximately_equal.hpp"
#include "algebra/impl/array_soa_boolean.hpp"
#include "algebra/impl/array_soa_casts.hpp"
#include "algebra/impl/array_soa_constants.hpp"
#include "algebra/impl/array_soa_math.hpp"
#include "algebra/impl/array_soa_types.hpp"
#include "detray/algebra/common/impl/matrix_getter.hpp"
// #include "detray/algebra/common/common.hpp"
#include "detray/algebra/generic/impl/generic_matrix.hpp"
#include "detray/algebra/generic/impl/generic_vector.hpp"

// System include(s).
#include <bit>
#include <cstddef>

namespace detray {

// The operators on @c algebra::storage::vector types are found by argument
// dependent lookup. They are not pulled into the detray namespace, because they
// would be ambiguous with the operators of the AoS array plugin

/// Define the plugin types
///
/// @tparam V Value type (scalar precision) of a lane.
/// @tparam W Number of SIMD lanes, must be a power of two.
/// @{
template <concepts::value V, std::size_t W = 8u>
struct array_soa {
  static_assert(std::has_single_bit(W),
                "The number of SIMD lanes must be a power of two");

  /// @returns the number of SIMD lanes
  static consteval std::size_t size() { return W; }

  /// Define scalar precision
  using value_type = V;

  template <concepts::element T>
  using simd = algebra::array_soa::simd<T, W>;

  /// Lane mask: one integer per lane with the same size as the values
  using boolean = typename simd<value_type>::mask_type;

  /// Linear Algebra type definitions
  /// @{
  using scalar = simd<value_type>;
  using index_type = algebra::array_soa::index_type;
  using transform3D = algebra::array_soa::transform3<value_type, W>;
  using point2D = algebra::array_soa::point2<value_type, W>;
  using point3D = algebra::array_soa::point3<value_type, W>;
  using vector2D = algebra::array_soa::vector2<value_type, W>;
  using vector3D = algebra::array_soa::vector3<value_type, W>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::array_soa::matrix_type<value_type, ROWS, COLS, W>;
  /// @}
};
/// @}

namespace getter {

/// @name Getter functions on @c algebra::array_soa types
/// @{

using algebra::common::getter::block;
using algebra::common::getter::element;
using algebra::common::getter::set_block;
using algebra::common::getter::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::array_soa types
/// @{

using algebra::generic::math::cross;
using algebra::generic::math::dot;
using algebra::generic::math::eta;
using algebra::generic::math::norm;
using algebra::generic::math::normalize;
using algebra::generic::math::perp;
using algebra::generic::math::phi;
using algebra::generic::math::theta;

/// @}

}  // namespace vector

// Produces clash with matrix typedefs in other plugins
namespace matrix {

using algebra::common::determinant;
using algebra::common::identity;
using algebra::common::inverse;
using algebra::common::set_identity;
using algebra::common::set_zero;
using algebra::common::transpose;
using algebra::common::zero;

using algebra::generic::math::cholesky_decomposition;
using algebra::generic::math::column_wise_cross;
using algebra::generic::math::column_wise_multiply;
using algebra::generic::math::cross_matrix;
using algebra::generic::math::outer_product;
using algebra::generic::math::set_inplace_product_left;
using algebra::generic::math::set_inplace_product_left_transpose;
using algebra::generic::math::set_inplace_product_right;
using algebra::generic::math::set_inplace_product_right_transpose;
using algebra::generic::math::set_product;
using algebra::generic::math::set_product_left_transpose;
using algebra::generic::math::set_product_right_transpose;
using algebra::generic::math::transposed_product;

}  // namespace matrix

}  // namespace detray
