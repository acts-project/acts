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
#include "algebra/impl/vc_aos_getter.hpp"
#include "algebra/impl/vc_aos_matrix.hpp"
#include "algebra/impl/vc_aos_types.hpp"
#include "algebra/impl/vc_aos_vector.hpp"
#include "detray/algebra/common/boolean.hpp"
#include "detray/algebra/generic/generic.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

namespace detray {

/// Define the plugin types
/// @{
template <concepts::value V>
struct vc_aos {
  /// Define scalar type
  using value_type = V;

  template <concepts::element T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using index_type = algebra::vc_aos::index_type;
  using transform3D = algebra::vc_aos::transform3<value_type>;
  using point2D = algebra::vc_aos::point2<value_type>;
  using point3D = algebra::vc_aos::point3<value_type>;
  using vector2D = algebra::vc_aos::vector2<value_type>;
  using vector3D = algebra::vc_aos::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::vc_aos::matrix_type<value_type, ROWS, COLS>;
};
/// @}

namespace getter {

/// @name Getter functions on @c algebra::vc_aos::matrix_type
/// @{

using algebra::vc_aos::storage::block;
using algebra::vc_aos::storage::element;
using algebra::vc_aos::storage::set_block;
using algebra::vc_aos::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc_aos types
/// @{

// Vc array specific
using algebra::vc_aos::math::cross;
using algebra::vc_aos::math::dot;
using algebra::vc_aos::math::eta;
using algebra::vc_aos::math::norm;
using algebra::vc_aos::math::normalize;
using algebra::vc_aos::math::perp;
using algebra::vc_aos::math::phi;
using algebra::vc_aos::math::theta;

/// @}

}  // namespace vector

// Use special algorithms for 4 dimensional matrices
namespace algebra::generic {

// Determinant algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct determinant_selector<4, vc_aos::matrix_type<T, ROWS, COLS>,
                            vc_aos::element_getter> {
  using type =
      matrix::determinant::hard_coded<vc_aos::matrix_type<T, ROWS, COLS>>;
};

// Inversion algorithms
template <concepts::scalar T, auto ROWS, auto COLS>
struct inversion_selector<4, vc_aos::matrix_type<T, ROWS, COLS>,
                          vc_aos::element_getter> {
  using type = matrix::inverse::hard_coded<vc_aos::matrix_type<T, ROWS, COLS>>;
};

}  // namespace algebra::generic

namespace matrix {

/// @name Matrix functions on @c algebra::vc_aos types
/// @{

using algebra::vc_aos::math::determinant;
using algebra::vc_aos::math::identity;
using algebra::vc_aos::math::inverse;
using algebra::vc_aos::math::set_identity;
using algebra::vc_aos::math::set_zero;
using algebra::vc_aos::math::transpose;
using algebra::vc_aos::math::zero;

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
