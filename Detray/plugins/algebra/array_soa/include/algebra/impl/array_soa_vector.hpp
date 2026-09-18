// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/array_soa_math.hpp"
#include "detray/algebra/common/vector.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::algebra::array_soa::math {

/// This method retrieves phi from a vector, vector base with rows >= 2
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam W number of lanes
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t, std::size_t W,
          template <typename, std::size_t> class array_t>
  requires(N >= 2)
DETRAY_HOST_DEVICE constexpr auto phi(
    const algebra::storage::vector<N, array_soa::simd<value_t, W>, array_t>
        &v) {
  return algebra::math::atan2(v[1], v[0]);
}

/// This method retrieves the perpendicular magnitude of a vector with rows >= 2
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam W number of lanes
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t, std::size_t W,
          template <typename, std::size_t> class array_t>
  requires(N >= 2)
DETRAY_HOST_DEVICE constexpr auto perp(
    const algebra::storage::vector<N, array_soa::simd<value_t, W>, array_t>
        &v) {
  return algebra::math::sqrt(algebra::math::fma(v[0], v[0], v[1] * v[1]));
}

/// This method retrieves theta from a vector, vector base with rows >= 3
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam W number of lanes
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t, std::size_t W,
          template <typename, std::size_t> class array_t>
  requires(N >= 3)
DETRAY_HOST_DEVICE constexpr auto theta(
    const algebra::storage::vector<N, array_soa::simd<value_t, W>, array_t>
        &v) {
  return algebra::math::atan2(perp(v), v[2]);
}

/// Cross product between two input vectors - 3 Dim
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam W number of lanes
/// @tparam array_t array type that holds the vector elements
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector (expression) representing the cross product
template <std::size_t N, concepts::value value_t, std::size_t W,
          template <typename, std::size_t> class array_t>
  requires(N == 3)
DETRAY_HOST_DEVICE constexpr algebra::storage::vector<
    N, array_soa::simd<value_t, W>, array_t>
cross(
    const algebra::storage::vector<N, array_soa::simd<value_t, W>, array_t> &a,
    const algebra::storage::vector<N, array_soa::simd<value_t, W>, array_t>
        &b) {
  return {algebra::math::fma(a[1], b[2], -b[1] * a[2]),
          algebra::math::fma(a[2], b[0], -b[2] * a[0]),
          algebra::math::fma(a[0], b[1], -b[0] * a[1])};
}

/// Dot product between two input vectors
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam W number of lanes
/// @tparam array_t array type that holds the vector elements
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <std::size_t N, concepts::value value_t, std::size_t W,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr array_soa::simd<value_t, W> dot(
    const algebra::storage::vector<N, array_soa::simd<value_t, W>, array_t> &a,
    const algebra::storage::vector<N, array_soa::simd<value_t, W>, array_t>
        &b) {
  auto ret = a[0] * b[0];

  for (unsigned int i{1u}; i < N; i++) {
    ret = algebra::math::fma(a[i], b[i], ret);
  }

  return ret;
}

/// This method retrieves the norm of a vector, no dimension restriction
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam W number of lanes
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t, std::size_t W,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr auto norm(
    const algebra::storage::vector<N, array_soa::simd<value_t, W>, array_t>
        &v) {
  return algebra::math::sqrt(dot(v, v));
}

/// Get a normalized version of the input vector
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam W number of lanes
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t, std::size_t W,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::vector<
    N, array_soa::simd<value_t, W>, array_t>
normalize(const algebra::storage::vector<N, array_soa::simd<value_t, W>,
                                         array_t> &v) {
  return (array_soa::simd<value_t, W>::One() / norm(v)) * v;
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @tparam N dimension of the vector
/// @tparam value_t value type in the simd vectors
/// @tparam W number of lanes
/// @tparam array_t array type that holds the vector elements
///
/// @param v the input vector
template <std::size_t N, concepts::value value_t, std::size_t W,
          template <typename, std::size_t> class array_t>
  requires(N >= 3)
DETRAY_HOST_DEVICE constexpr auto eta(
    const algebra::storage::vector<N, array_soa::simd<value_t, W>, array_t>
        &v) {
  return algebra::math::atanh(v[2] / norm(v));
}

}  // namespace detray::algebra::array_soa::math
