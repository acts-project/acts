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
#include "algebra/impl/vc_soa_getter.hpp"
#include "algebra/impl/vc_soa_math.hpp"
#include "algebra/impl/vc_soa_matrix.hpp"
#include "algebra/impl/vc_soa_types.hpp"
#include "algebra/impl/vc_soa_vector.hpp"
#include "detray/algebra/generic/impl/generic_matrix.hpp"

// System include(s).
#include <cassert>
#include <type_traits>

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
template <concepts::value V>
struct vc_soa {
  /// Define scalar precision
  using value_type = V;

  template <concepts::element T>
  using simd = Vc::Vector<T>;

  using boolean = Vc::Mask<V>;

  /// Linear Algebra type definitions
  /// @{
  using scalar = simd<value_type>;
  using index_type = algebra::vc_soa::index_type;
  using transform3D = algebra::vc_soa::transform3<value_type>;
  using point2D = algebra::vc_soa::point2<value_type>;
  using point3D = algebra::vc_soa::point3<value_type>;
  using vector2D = algebra::vc_soa::vector2<value_type>;
  using vector3D = algebra::vc_soa::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::vc_soa::matrix_type<value_type, ROWS, COLS>;
  /// @}
};
/// @}

namespace getter {

/// @name Getter functions on @c algebra::vc_soa types
/// @{

using algebra::vc_soa::storage::block;
using algebra::vc_soa::storage::element;
using algebra::vc_soa::storage::set_block;
using algebra::vc_soa::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::vc_soa types
/// @{

using algebra::vc_soa::math::cross;
using algebra::vc_soa::math::dot;
using algebra::vc_soa::math::eta;
using algebra::vc_soa::math::norm;
using algebra::vc_soa::math::normalize;
using algebra::vc_soa::math::perp;
using algebra::vc_soa::math::phi;
using algebra::vc_soa::math::theta;

/// @}

}  // namespace vector

// Produces clash with matrix typedefs in other plugins
namespace matrix {

using algebra::vc_soa::math::determinant;
using algebra::vc_soa::math::identity;
using algebra::vc_soa::math::inverse;
using algebra::vc_soa::math::set_identity;
using algebra::vc_soa::math::set_zero;
using algebra::vc_soa::math::transpose;
using algebra::vc_soa::math::zero;

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
