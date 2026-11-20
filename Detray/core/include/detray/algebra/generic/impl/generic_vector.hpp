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

// The math functions are called without namespace qualifier, after a
// using-declaration of the @c algebra::math overload set. With this, argument
// dependent lookup finds the overloads for SIMD scalar types at the point of
// instantiation. This header can then be included before the math overloads of
// a plugin are declared (e.g. when several plugins are used together).
//
// Note: A plugin with a SIMD scalar type has to make its math overloads
// visible in the namespace of the SIMD type.

namespace detray::algebra::generic::math {

/// This method retrieves phi from a vector with rows >= 2
///
/// @param v the input vector
template <concepts::vector vector_t>
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector_t> phi(
    const vector_t &v) noexcept {
  using element_getter_t = detray::traits::element_getter_t<vector_t>;
  using algebra::math::atan2;

  return atan2(element_getter_t{}(v, 1), element_getter_t{}(v, 0));
}

/// This method retrieves the perpendicular magnitude of a vector with rows >= 2
///
/// @param v the input vector
template <concepts::vector vector_t>
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector_t> perp(
    const vector_t &v) noexcept {
  using element_getter_t = detray::traits::element_getter_t<vector_t>;
  using algebra::math::fma;
  using algebra::math::sqrt;

  return sqrt(fma(element_getter_t{}(v, 0), element_getter_t{}(v, 0),
                  element_getter_t{}(v, 1) * element_getter_t{}(v, 1)));
}

/// This method retrieves theta from a vector with rows >= 3
///
/// @param v the input vector
template <concepts::vector vector_t>
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector_t> theta(
    const vector_t &v) noexcept {
  using element_getter_t = detray::traits::element_getter_t<vector_t>;
  using algebra::math::atan2;

  return atan2(perp(v), element_getter_t{}(v, 2));
}

/// Cross product between two input vectors - 3 Dim
///
/// @tparam vector1_t first vector or column matrix type
/// @tparam vector2_t second vector or column matrix type
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector representing the cross product
template <typename vector1_t, typename vector2_t>
  requires(
      (concepts::vector3D<vector1_t> || concepts::column_matrix3D<vector1_t>) &&
      (concepts::vector3D<vector2_t> || concepts::column_matrix3D<vector2_t>))
DETRAY_HOST_DEVICE constexpr detray::traits::vector_t<vector1_t> cross(
    const vector1_t &a, const vector2_t &b) {
  using element_getter_t = detray::traits::element_getter_t<vector1_t>;
  using algebra::math::fma;

  return {fma(element_getter_t{}(a, 1), element_getter_t{}(b, 2),
              -element_getter_t{}(b, 1) * element_getter_t{}(a, 2)),
          fma(element_getter_t{}(a, 2), element_getter_t{}(b, 0),
              -element_getter_t{}(b, 2) * element_getter_t{}(a, 0)),
          fma(element_getter_t{}(a, 0), element_getter_t{}(b, 1),
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
/// @return the scalar dot product value
template <typename vector1_t, typename vector2_t>
  requires((concepts::vector<vector1_t> ||
            concepts::column_matrix<vector1_t>) &&
           (concepts::vector<vector2_t> || concepts::column_matrix<vector2_t>))
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector1_t> dot(
    const vector1_t &a, const vector2_t &b) {
  using scalar_t = detray::traits::scalar_t<vector1_t>;
  using index_t = detray::traits::index_t<vector1_t>;
  using element_getter_t = detray::traits::element_getter_t<vector1_t>;
  using algebra::math::fma;

  scalar_t ret{element_getter_t{}(a, 0) * element_getter_t{}(b, 0)};

  for (index_t i = 1; i < detray::traits::size<vector1_t>; i++) {
    ret = fma(element_getter_t{}(a, i), element_getter_t{}(b, i), ret);
  }

  return ret;
}

/// This method retrieves the norm of a vector with rows >= 3
///
/// @param v the input vector
template <concepts::vector vector_t>
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector_t> norm(
    const vector_t &v) {
  using algebra::math::sqrt;

  return sqrt(dot(v, v));
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @param v the input vector
template <concepts::vector vector_t>
DETRAY_HOST_DEVICE constexpr detray::traits::scalar_t<vector_t> eta(
    const vector_t &v) noexcept {
  using element_getter_t = detray::traits::element_getter_t<vector_t>;
  using algebra::math::atanh;

  return atanh(element_getter_t{}(v, 2) / norm(v));
}

/// Get a normalized version of the input vector
///
/// @tparam vector_t vector or column matrix type
///
/// @param v the input vector
///
/// @returns the normalized vector
template <concepts::vector vector_t>
DETRAY_HOST_DEVICE constexpr vector_t normalize(const vector_t &v) {
  using scalar_t = detray::traits::scalar_t<vector_t>;

  // One division and N multiplications instead of N divisions
  return (static_cast<scalar_t>(1.f) / norm(v)) * v;
}

}  // namespace detray::algebra::generic::math
