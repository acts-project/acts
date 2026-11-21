// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/type_traits.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// Generic vector implementations
//
// A plugin opts in to these functions for its vector (and column matrix) types:
// it specializes @c traits::plugin_uses_generic_vector. The functions only
// match the opted-in types.
//
// Note: A plugin with a SIMD scalar type has to declare its math overloads in
// @c algebra::math and include them before this header (see
// detray/definitions/algebra.hpp).

namespace detray::algebra::generic::math {

/// @brief Retrieves phi from a vector with rows >= 2
///
/// @tparam vector_t vector type that uses the generic implementation
///
/// @param v the input vector
///
/// @returns the azimuthal angle phi
template <concepts::generic_vector vector_t>
  requires(detray::traits::size<vector_t> >= 2)
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector_t> phi(
    const vector_t &v) noexcept {
  using element_getter_t = detray::traits::element_getter_t<vector_t>;

  return algebra::math::atan2(element_getter_t{}(v, 1),
                              element_getter_t{}(v, 0));
}

/// @brief Retrieves the perpendicular magnitude of a vector with rows >= 2
///
/// @tparam vector_t vector type that uses the generic implementation
///
/// @param v the input vector
///
/// @returns the transverse component perp
template <concepts::generic_vector vector_t>
  requires(detray::traits::size<vector_t> >= 2)
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector_t> perp(
    const vector_t &v) noexcept {
  using element_getter_t = detray::traits::element_getter_t<vector_t>;

  return algebra::math::sqrt(
      algebra::math::fma(element_getter_t{}(v, 0), element_getter_t{}(v, 0),
                         element_getter_t{}(v, 1) * element_getter_t{}(v, 1)));
}

/// @brief Retrieves theta from a vector with rows >= 3
///
/// @tparam vector_t vector type that uses the generic implementation
///
/// @param v the input vector
///
/// @returns the polar angle theta
template <concepts::generic_vector vector_t>
  requires(detray::traits::size<vector_t> >= 3)
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector_t> theta(
    const vector_t &v) noexcept {
  using element_getter_t = detray::traits::element_getter_t<vector_t>;

  return algebra::math::atan2(algebra::generic::math::perp(v),
                              element_getter_t{}(v, 2));
}

/// Cross product between two input vectors - 3 Dim
///
/// @tparam vector1_t first vector or column matrix type
/// @tparam vector2_t second vector or column matrix type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @returns the cross product of @p a and @p b
template <typename vector1_t, typename vector2_t>
  requires(
      (concepts::vector3D<vector1_t> || concepts::column_matrix3D<vector1_t>) &&
      (concepts::vector3D<vector2_t> || concepts::column_matrix3D<vector2_t>) &&
      detray::traits::plugin_uses_generic_vector_v<vector1_t> &&
      detray::traits::plugin_uses_generic_vector_v<vector2_t>)
DETRAY_HOST_DEVICE constexpr detray::traits::vector_t<vector1_t> cross(
    const vector1_t &a, const vector2_t &b) {
  using element_getter_t = detray::traits::element_getter_t<vector1_t>;

  return {
      algebra::math::fma(element_getter_t{}(a, 1), element_getter_t{}(b, 2),
                         -element_getter_t{}(b, 1) * element_getter_t{}(a, 2)),
      algebra::math::fma(element_getter_t{}(a, 2), element_getter_t{}(b, 0),
                         -element_getter_t{}(b, 2) * element_getter_t{}(a, 0)),
      algebra::math::fma(element_getter_t{}(a, 0), element_getter_t{}(b, 1),
                         -element_getter_t{}(b, 0) * element_getter_t{}(a, 1))};
}

/// Dot product between two input vectors
///
/// @tparam vector1_t first vector or column matrix type
/// @tparam vector2_t second vector or column matrix type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @returns the dot product of @p a and @p b
template <typename vector1_t, typename vector2_t>
  requires(
      (concepts::vector<vector1_t> || concepts::column_matrix<vector1_t>) &&
      (concepts::vector<vector2_t> || concepts::column_matrix<vector2_t>) &&
      detray::traits::plugin_uses_generic_vector_v<vector1_t> &&
      detray::traits::plugin_uses_generic_vector_v<vector2_t>)
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector1_t> dot(
    const vector1_t &a, const vector2_t &b) {
  using scalar_t = detray::traits::scalar_t<vector1_t>;
  using index_t = detray::traits::index_t<vector1_t>;
  using element_getter_t = detray::traits::element_getter_t<vector1_t>;

  scalar_t ret{element_getter_t{}(a, 0) * element_getter_t{}(b, 0)};

  for (index_t i = 1; i < detray::traits::size<vector1_t>; i++) {
    ret = algebra::math::fma(element_getter_t{}(a, i), element_getter_t{}(b, i),
                             ret);
  }

  return ret;
}

/// @brief Retrieves the norm of a vector
///
/// @tparam vector_t vector type that uses the generic implementation
///
/// @param v the input vector
///
/// @returns the euclidean norm of @p v
template <concepts::generic_vector vector_t>
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector_t> norm(
    const vector_t &v) {
  return algebra::math::sqrt(algebra::generic::math::dot(v, v));
}

/// @brief Retrieves the pseudo-rapidity from a vector
///
/// @tparam vector_t vector type that uses the generic implementation
///
/// @param v the input vector
///
/// @returns the pseudorapidity eta
template <concepts::generic_vector vector_t>
  requires(detray::traits::size<vector_t> >= 3)
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector_t> eta(
    const vector_t &v) noexcept {
  using element_getter_t = detray::traits::element_getter_t<vector_t>;

  return algebra::math::atanh(element_getter_t{}(v, 2) /
                              algebra::generic::math::norm(v));
}

/// @brief Retrieves a normalized version of the input vector
///
/// @tparam vector_t vector type that uses the generic implementation
///
/// @param v the input vector
///
/// @returns the unit vector in the direction of @p v
template <concepts::generic_vector vector_t>
DETRAY_HOST_DEVICE constexpr vector_t normalize(const vector_t &v) {
  using scalar_t = detray::traits::scalar_t<vector_t>;

  // One division and N multiplications instead of N divisions
  return (static_cast<scalar_t>(1.f) / algebra::generic::math::norm(v)) * v;
}

}  // namespace detray::algebra::generic::math
