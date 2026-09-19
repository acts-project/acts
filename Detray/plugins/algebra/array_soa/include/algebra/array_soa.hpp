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
#include "algebra/impl/array_soa_getter.hpp"
#include "algebra/impl/array_soa_math.hpp"
#include "algebra/impl/array_soa_matrix.hpp"
#include "algebra/impl/array_soa_types.hpp"
#include "algebra/impl/array_soa_vector.hpp"
#include "detray/algebra/generic/impl/generic_matrix.hpp"

// System include(s).
#include <cstddef>

namespace detray {

/// @name Operators on @c algebra::storage::vector types
/// @{

using algebra::storage::operator*;
using algebra::storage::operator/;
using algebra::storage::operator-;
using algebra::storage::operator+;

/// @}

/// Define the plugin types
/// @{
template <concepts::value V, std::size_t W = 8u>
struct array_soa {
  /// @returns the widh of the SIMD lane
  static consteval std::size_t size() { return W; }

  /// Define scalar precision
  using value_type = V;

  template <concepts::element T>
  using simd = algebra::array_soa::simd<T, W>;

  using boolean = algebra::array_soa::mask<W>;

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

using algebra::array_soa::storage::block;
using algebra::array_soa::storage::element;
using algebra::array_soa::storage::set_block;
using algebra::array_soa::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::array_soa types
/// @{

using algebra::array_soa::math::cross;
using algebra::array_soa::math::dot;
using algebra::array_soa::math::eta;
using algebra::array_soa::math::norm;
using algebra::array_soa::math::normalize;
using algebra::array_soa::math::perp;
using algebra::array_soa::math::phi;
using algebra::array_soa::math::theta;

/// @}

}  // namespace vector

// Produces clash with matrix typedefs in other plugins
namespace matrix {

using algebra::array_soa::math::determinant;
using algebra::array_soa::math::identity;
using algebra::array_soa::math::inverse;
using algebra::array_soa::math::set_identity;
using algebra::array_soa::math::set_zero;
using algebra::array_soa::math::transpose;
using algebra::array_soa::math::zero;

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
