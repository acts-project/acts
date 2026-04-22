// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/impl/eigen_getter.hpp"
#include "algebra/impl/eigen_matrix.hpp"
#include "algebra/impl/eigen_types.hpp"
#include "algebra/impl/eigen_vector.hpp"
#include "detray/algebra/common/boolean.hpp"
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/generic/generic.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC

namespace detray {

/// Define the plugin types
/// @{
template <concepts::value V>
struct eigen {
  /// Define scalar type
  using value_type = V;

  template <concepts::element T>
  using simd = T;

  using boolean = bool;
  using scalar = value_type;
  using index_type = algebra::eigen::index_type;
  using transform3D = algebra::eigen::transform3<value_type>;
  using point2D = algebra::eigen::point2<value_type>;
  using point3D = algebra::eigen::point3<value_type>;
  using vector2D = algebra::eigen::vector2<value_type>;
  using vector3D = algebra::eigen::vector3<value_type>;

  template <std::size_t ROWS, std::size_t COLS>
  using matrix = algebra::eigen::matrix_type<value_type, ROWS, COLS>;
};
/// @}

namespace getter {

/// @name Getter functions on @c algebra::eigen
/// @{

using algebra::eigen::storage::block;
using algebra::eigen::storage::element;
using algebra::eigen::storage::set_block;
using algebra::eigen::storage::vector;

/// @}

}  // namespace getter

namespace vector {

/// @name Vector functions on @c algebra::eigen::storage_type
/// @{

using algebra::eigen::math::cross;
using algebra::eigen::math::dot;
using algebra::eigen::math::eta;
using algebra::eigen::math::norm;
using algebra::eigen::math::normalize;
using algebra::eigen::math::perp;
using algebra::eigen::math::phi;
using algebra::eigen::math::theta;

}  // namespace vector

namespace matrix {

/// @name Matrix functions on @c algebra::eigen::storage_type
/// @{

using algebra::eigen::math::cholesky_decomposition;
using algebra::eigen::math::column_wise_cross;
using algebra::eigen::math::column_wise_multiply;
using algebra::eigen::math::determinant;
using algebra::eigen::math::identity;
using algebra::eigen::math::inverse;
using algebra::eigen::math::outer_product;
using algebra::eigen::math::set_identity;
using algebra::eigen::math::set_zero;
using algebra::eigen::math::transpose;
using algebra::eigen::math::zero;

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
