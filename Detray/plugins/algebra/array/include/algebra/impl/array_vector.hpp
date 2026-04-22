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
#include "algebra/impl/array_operators.hpp"
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/generic/generic.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::array {

/// This method retrieves phi from a vector with rows >= 2
///
/// @param v the input vector
template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
  requires(N >= 2)
DETRAY_HOST_DEVICE constexpr scalar_t phi(const array_t<scalar_t, N> &v) {
  return algebra::generic::math::phi(v);
}

/// This method retrieves the perpendicular magnitude of a vector with rows >= 2
///
/// @param v the input vector
template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
  requires(N >= 2)
DETRAY_HOST_DEVICE constexpr scalar_t perp(const array_t<scalar_t, N> &v) {
  return algebra::generic::math::perp(v);
}

/// This method retrieves theta from a vector with rows >= 3
///
/// @param v the input vector
template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
  requires(N >= 2)
DETRAY_HOST_DEVICE constexpr scalar_t theta(const array_t<scalar_t, N> &v) {
  return algebra::generic::math::theta(v);
}

/// Cross product between two input vectors - 3 Dim
///
/// @tparam index_t the index type for this plugin
/// @tparam array_t the array type the plugin is based on
/// @tparam scalar_t the scalar type
/// @tparam N the dimension of the vectors (minimum 3)
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector representing the cross product
/// @{
template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
  requires(N >= 3)
DETRAY_HOST_DEVICE constexpr array_t<scalar_t, N> cross(
    const array_t<scalar_t, N> &a, const array_t<scalar_t, N> &b) {
  return algebra::generic::math::cross(a, b);
}

template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
  requires(N >= 3)
DETRAY_HOST_DEVICE constexpr array_t<scalar_t, N> cross(
    const array_t<scalar_t, N> &a, const array_t<array_t<scalar_t, N>, 1> &b) {
  return algebra::generic::math::cross(a, b);
}

template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
  requires(N >= 3)
DETRAY_HOST_DEVICE constexpr array_t<scalar_t, N> cross(
    const array_t<array_t<scalar_t, N>, 1> &a, const array_t<scalar_t, N> &b) {
  return algebra::generic::math::cross(a, b);
}

template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
  requires(N >= 3)
DETRAY_HOST_DEVICE constexpr array_t<scalar_t, N> cross(
    const array_t<array_t<scalar_t, N>, 1> &a,
    const array_t<array_t<scalar_t, N>, 1> &b) {
  return algebra::generic::math::cross(a, b);
}
/// @}

/// Dot product between two input vectors/column matrices
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
/// @{
template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
DETRAY_HOST_DEVICE constexpr scalar_t dot(const array_t<scalar_t, N> &a,
                                          const array_t<scalar_t, N> &b) {
  array_t<scalar_t, N> tmp;
  for (index_t i = 0; i < N; i++) {
    tmp[i] = a[i] * b[i];
  }
  scalar_t ret{0.f};
  for (index_t i = 0; i < N; i++) {
    ret += tmp[i];
  }
  return ret;
}

template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
DETRAY_HOST_DEVICE constexpr scalar_t dot(
    const array_t<scalar_t, N> &a, const array_t<array_t<scalar_t, N>, 1> &b) {
  array_t<scalar_t, N> tmp;
  for (index_t i = 0; i < N; i++) {
    tmp[i] = a[i] * b[0][i];
  }
  scalar_t ret{0.f};
  for (index_t i = 0; i < N; i++) {
    ret += tmp[i];
  }
  return ret;
}

template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
DETRAY_HOST_DEVICE constexpr scalar_t dot(
    const array_t<array_t<scalar_t, N>, 1> &a, const array_t<scalar_t, N> &b) {
  array_t<scalar_t, N> tmp;
  for (index_t i = 0; i < N; i++) {
    tmp[i] = a[0][i] * b[i];
  }
  scalar_t ret{0.f};
  for (index_t i = 0; i < N; i++) {
    ret += tmp[i];
  }
  return ret;
}

template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
DETRAY_HOST_DEVICE constexpr scalar_t dot(
    const array_t<array_t<scalar_t, N>, 1> &a,
    const array_t<array_t<scalar_t, N>, 1> &b) {
  array_t<scalar_t, N> tmp;
  for (index_t i = 0; i < N; i++) {
    tmp[i] = a[0][i] * b[0][i];
  }
  scalar_t ret{0.f};
  for (index_t i = 0; i < N; i++) {
    ret += tmp[i];
  }
  return ret;
}
/// @}

/// This method retrieves the norm of a vector, no dimension restriction
///
/// @param v the input vector
template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
  requires(N >= 2)
DETRAY_HOST_DEVICE constexpr scalar_t norm(const array_t<scalar_t, N> &v) {
  return algebra::math::sqrt(dot(v, v));
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @param v the input vector
template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
  requires(N >= 3)
DETRAY_HOST_DEVICE constexpr scalar_t eta(
    const array_t<scalar_t, N> &v) noexcept {
  return algebra::math::atanh(v[2] / norm(v));
}

/// Get a normalized version of the input vector
///
/// @param v the input vector
template <concepts::index index_t, template <typename, index_t> class array_t,
          concepts::scalar scalar_t, index_t N>
DETRAY_HOST_DEVICE constexpr array_t<scalar_t, N> normalize(
    const array_t<scalar_t, N> &v) {
  return (static_cast<scalar_t>(1.) / norm(v)) * v;
}

}  // namespace detray::algebra::array
