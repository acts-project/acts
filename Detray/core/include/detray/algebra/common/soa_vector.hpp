// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/common/constants.hpp"
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/common/vector.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <cstddef>

/// Vector functions for the SoA plugins (Vc SoA and array SoA)
///
/// The functions operate on @c algebra::storage::vector types with a SIMD
/// scalar. All lanes are evaluated at the same time.
///
/// @note The plugin has to declare the overloads of @c algebra::math::fma,
/// @c sqrt, @c atan2, @c atanh and of @c algebra::constants::one for its SIMD
/// type before this header is included.
namespace detray::algebra::soa {

/// Get the azimuthal angle of a vector in every lane
///
/// @tparam N Number of vector elements
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @param v the input vector
///
/// @returns the azimuthal angle phi
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
  requires(N >= 2)
DETRAY_HOST_DEVICE constexpr scalar_t phi(
    const algebra::storage::vector<N, scalar_t, array_t> &v) {
  return algebra::math::atan2(v[1], v[0]);
}

/// Get the transverse component of a vector in every lane
///
/// @tparam N Number of vector elements
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @param v the input vector
///
/// @returns the transverse component perp
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
  requires(N >= 2)
DETRAY_HOST_DEVICE constexpr scalar_t perp(
    const algebra::storage::vector<N, scalar_t, array_t> &v) {
  return algebra::math::sqrt(algebra::math::fma(v[0], v[0], v[1] * v[1]));
}

/// Get the polar angle of a vector in every lane
///
/// @tparam N Number of vector elements
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @param v the input vector
///
/// @returns the polar angle theta
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
  requires(N >= 3)
DETRAY_HOST_DEVICE constexpr scalar_t theta(
    const algebra::storage::vector<N, scalar_t, array_t> &v) {
  return algebra::math::atan2(perp(v), v[2]);
}

/// Get the cross product of two 3D vectors in every lane
///
/// @tparam N Number of vector elements
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @returns the cross product of @param a and @param b
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
  requires(N == 3)
DETRAY_HOST_DEVICE constexpr algebra::storage::vector<N, scalar_t, array_t>
cross(const algebra::storage::vector<N, scalar_t, array_t> &a,
      const algebra::storage::vector<N, scalar_t, array_t> &b) {
  return {algebra::math::fma(a[1], b[2], -b[1] * a[2]),
          algebra::math::fma(a[2], b[0], -b[2] * a[0]),
          algebra::math::fma(a[0], b[1], -b[0] * a[1])};
}

/// Get the dot product of two vectors in every lane
///
/// @tparam N Number of vector elements
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @returns the dot product of @param a and @param b
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr scalar_t dot(
    const algebra::storage::vector<N, scalar_t, array_t> &a,
    const algebra::storage::vector<N, scalar_t, array_t> &b) {
  scalar_t ret{a[0] * b[0]};

  for (std::size_t i = 1u; i < N; ++i) {
    ret = algebra::math::fma(a[i], b[i], ret);
  }

  return ret;
}

/// Get the norm of a vector in every lane
///
/// @tparam N Number of vector elements
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @param v the input vector
///
/// @returns the euclidean norm of @param v
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr scalar_t norm(
    const algebra::storage::vector<N, scalar_t, array_t> &v) {
  return algebra::math::sqrt(dot(v, v));
}

/// Get the normalized vector in every lane
///
/// @tparam N Number of vector elements
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @param v the input vector
///
/// @returns the unit vector in the direction of @param v
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::vector<N, scalar_t, array_t>
normalize(const algebra::storage::vector<N, scalar_t, array_t> &v) {
  return (algebra::constants::one<scalar_t>() / norm(v)) * v;
}

/// Get the pseudorapidity of a vector in every lane
///
/// @tparam N Number of vector elements
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @param v the input vector
///
/// @returns the pseudorapidity eta
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
  requires(N >= 3)
DETRAY_HOST_DEVICE constexpr scalar_t eta(
    const algebra::storage::vector<N, scalar_t, array_t> &v) {
  return algebra::math::atanh(v[2] / norm(v));
}

}  // namespace detray::algebra::soa
