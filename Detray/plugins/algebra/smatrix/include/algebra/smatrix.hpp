// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/smatrix_getter.hpp"
#include "algebra/impl/smatrix_matrix.hpp"
#include "algebra/impl/smatrix_types.hpp"
#include "algebra/impl/smatrix_vector.hpp"
#include "detray/algebra/common/boolean.hpp"
#include "detray/algebra/generic/generic.hpp"

namespace detray {

/// Define the plugin types
/// @{
template <concepts::value V>
struct smatrix {
  /// Define scalar type
  using value_type = V;

  template <concepts::element T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using index_type = algebra::smatrix::index_type;
  using transform3D = algebra::smatrix::transform3<value_type>;
  using point2D = algebra::smatrix::point2<value_type>;
  using point3D = algebra::smatrix::point3<value_type>;
  using vector2D = algebra::smatrix::vector2<value_type>;
  using vector3D = algebra::smatrix::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::smatrix::matrix_type<value_type, ROWS, COLS>;
};
/// @}

namespace getter {

/// @name Getter functions on @c algebra::smatrix::matrix_type
/// @{

using algebra::smatrix::storage::block;
using algebra::smatrix::storage::element;
using algebra::smatrix::storage::set_block;
using algebra::smatrix::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::smatrix::storage_type
/// @{

using algebra::smatrix::math::cross;
using algebra::smatrix::math::dot;
using algebra::smatrix::math::eta;
using algebra::smatrix::math::norm;
using algebra::smatrix::math::normalize;
using algebra::smatrix::math::perp;
using algebra::smatrix::math::phi;
using algebra::smatrix::math::theta;

/// @}

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::smatrix::storage_type
/// @{

using algebra::smatrix::math::cholesky_decomposition;
using algebra::smatrix::math::column_wise_cross;
using algebra::smatrix::math::column_wise_multiply;
using algebra::smatrix::math::determinant;
using algebra::smatrix::math::identity;
using algebra::smatrix::math::inverse;
using algebra::smatrix::math::outer_product;
using algebra::smatrix::math::set_identity;
using algebra::smatrix::math::set_zero;
using algebra::smatrix::math::transpose;
using algebra::smatrix::math::zero;

using algebra::generic::math::cross_matrix;
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
