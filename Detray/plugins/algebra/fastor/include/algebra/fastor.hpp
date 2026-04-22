// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/fastor_getter.hpp"
#include "algebra/impl/fastor_matrix.hpp"
#include "algebra/impl/fastor_types.hpp"
#include "algebra/impl/fastor_vector.hpp"
#include "detray/algebra/common/boolean.hpp"
#include "detray/algebra/generic/generic.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

namespace detray {

template <std::size_t N, concepts::scalar S>
constexpr bool operator==(const Fastor::Tensor<S, N>& lhs,
                          const Fastor::Tensor<S, N>& rhs) {
  return Fastor::isequal(lhs, rhs, 0.f);
}

/// Define the plugin types
/// @{
template <concepts::value V>
struct fastor {
  /// Define scalar type
  using value_type = V;

  template <concepts::element T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using index_type = algebra::fastor::index_type;
  using transform3D = algebra::fastor::transform3<value_type>;
  using point2D = algebra::fastor::point2<value_type>;
  using point3D = algebra::fastor::point3<value_type>;
  using vector2D = algebra::fastor::vector2<value_type>;
  using vector3D = algebra::fastor::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::fastor::matrix_type<value_type, ROWS, COLS>;
};
/// @}

namespace getter {

/// @name Getter functions on @c algebra::fastor::matrix_type
/// @{

using algebra::fastor::storage::block;
using algebra::fastor::storage::element;
using algebra::fastor::storage::set_block;
using algebra::fastor::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::fastor::storage_type
/// @{

using algebra::fastor::math::cross;
using algebra::fastor::math::dot;
using algebra::fastor::math::eta;
using algebra::fastor::math::norm;
using algebra::fastor::math::normalize;
using algebra::fastor::math::perp;
using algebra::fastor::math::phi;
using algebra::fastor::math::theta;

/// @}

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::fastor::storage_type
/// @{

using algebra::fastor::math::column_wise_cross;
using algebra::fastor::math::column_wise_multiply;
using algebra::fastor::math::determinant;
using algebra::fastor::math::identity;
using algebra::fastor::math::inverse;
using algebra::fastor::math::outer_product;
using algebra::fastor::math::set_identity;
using algebra::fastor::math::set_zero;
using algebra::fastor::math::transpose;
using algebra::fastor::math::zero;

using algebra::generic::math::cholesky_decomposition;
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
