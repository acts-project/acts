// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/array_getter.hpp"
#include "algebra/impl/array_matrix.hpp"
#include "algebra/impl/array_operators.hpp"
#include "algebra/impl/array_types.hpp"
#include "algebra/impl/array_vector.hpp"
#include "detray/algebra/generic/generic.hpp"

namespace detray {

/// Define the plugin types
/// @{
template <concepts::value V>
struct array {
  /// Define scalar type
  using value_type = V;

  template <concepts::element T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using index_type = algebra::array::index_type;
  using transform3D = algebra::array::transform3<value_type>;
  using point2D = algebra::array::point2<value_type>;
  using point3D = algebra::array::point3<value_type>;
  using vector2D = algebra::array::vector2<value_type>;
  using vector3D = algebra::array::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::array::matrix_type<value_type, ROWS, COLS>;
};
/// @}

/// @name Operators on @c algebra::array::storage_type
/// @{

using algebra::array::operator*;
using algebra::array::operator-;
using algebra::array::operator+;

/// @}

namespace getter {

/// @name Getter functions on @c algebra::array::storage_type
/// @{

using algebra::array::storage::block;
using algebra::array::storage::element;
using algebra::array::storage::set_block;
using algebra::array::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::array::storage_type
/// @{

// array specific implementations
using algebra::array::dot;
using algebra::array::normalize;

// generic implementations
using algebra::array::cross;
using algebra::array::eta;
using algebra::array::norm;
using algebra::array::perp;
using algebra::array::phi;
using algebra::array::theta;

/// @}

}  // namespace vector

// Use special algorithms for 4 dimensional matrices
namespace algebra::generic {

// Determinant algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct determinant_selector<4, array::matrix_type<T, ROWS, COLS>,
                            array::element_getter> {
  using type =
      matrix::determinant::hard_coded<array::matrix_type<T, ROWS, COLS>>;
};

// Inversion algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct inversion_selector<4, array::matrix_type<T, ROWS, COLS>,
                          array::element_getter> {
  using type = matrix::inverse::hard_coded<array::matrix_type<T, ROWS, COLS>>;
};

}  // namespace algebra::generic

namespace matrix {

/// @name Matrix functions on @c algebra::array::storage_type
/// @{

using algebra::array::identity;
using algebra::array::set_identity;
using algebra::array::set_zero;
using algebra::array::zero;

// Uses generic implementation in the background
using algebra::array::determinant;
using algebra::array::inverse;
using algebra::array::transpose;

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

/// @}

}  // namespace matrix

}  // namespace detray
