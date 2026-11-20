// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/vc_aos_approximately_equal.hpp"
#include "algebra/impl/vc_soa_boolean.hpp"
#include "algebra/impl/vc_soa_casts.hpp"
#include "algebra/impl/vc_soa_constants.hpp"
#include "algebra/impl/vc_soa_math.hpp"
#include "algebra/impl/vc_soa_types.hpp"
#include "detray/algebra/common/matrix_getter.hpp"
#include "detray/algebra/common/soa_matrix.hpp"
#include "detray/algebra/common/soa_vector.hpp"
#include "detray/algebra/generic/impl/generic_matrix.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

namespace detray {

// The operators on @c algebra::storage::vector types are found by argument
// dependent lookup. They are not pulled into the detray namespace, because they
// would be ambiguous with the operators of the AoS array plugin

/// Define the plugin types
/// @{
template <concepts::value V>
struct vc_soa {
  /// Define scalar precision
  using value_type = V;

  /// Backend vectorized representation of an element type.
  template <concepts::element T>
  using simd = Vc::Vector<T>;

  /// Boolean type used for element-wise predicates.
  using boolean = Vc::Mask<V>;

  /// Linear Algebra type definitions
  /// @{
  using scalar = simd<value_type>;
  /// Index type used for matrix and vector access.
  using index_type = algebra::vc_soa::index_type;
  /// Three-dimensional transform type.
  using transform3D = algebra::vc_soa::transform3<value_type>;
  /// Two-dimensional point type.
  using point2D = algebra::vc_soa::point2<value_type>;
  /// Three-dimensional point type.
  using point3D = algebra::vc_soa::point3<value_type>;
  /// Two-dimensional vector type.
  using vector2D = algebra::vc_soa::vector2<value_type>;
  /// Three-dimensional vector type.
  using vector3D = algebra::vc_soa::vector3<value_type>;

  /// Fixed-size matrix type with the requested dimensions.
  template <std::size_t ROWS, std::size_t COLS>
  /// Fixed-size matrix type with the requested row and column counts.
  using matrix = algebra::vc_soa::matrix_type<value_type, ROWS, COLS>;
  /// @}

  /// @returns the number of SIMD lanes
  static consteval std::size_t size() { return simd<V>::Size; }
};
/// @}

namespace getter {

/// @name Getter functions on @c algebra::vc_soa types
/// @{

using algebra::storage::getter::block;
using algebra::storage::getter::element;
using algebra::storage::getter::set_block;
using algebra::storage::getter::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc_soa types
/// @{

using algebra::soa::cross;
using algebra::soa::dot;
using algebra::soa::eta;
using algebra::soa::norm;
using algebra::soa::normalize;
using algebra::soa::perp;
using algebra::soa::phi;
using algebra::soa::theta;

/// @}

}  // namespace vector

// Produces clash with matrix typedefs in other plugins
namespace matrix {

using algebra::soa::determinant;
using algebra::soa::identity;
using algebra::soa::inverse;
using algebra::soa::set_identity;
using algebra::soa::set_zero;
using algebra::soa::transpose;
using algebra::soa::zero;

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
