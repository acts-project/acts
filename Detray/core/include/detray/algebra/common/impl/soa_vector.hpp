// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/common/impl/storage_vector.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/generic/impl/generic_vector.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <cstddef>

// Vector functions for the SoA plugins (Vc SoA and array SoA)
//
// The functions operate on algebra::storage::vector types with a SIMD
// scalar. All lanes are evaluated at the same time. The implementation is in
// algebra::generic::math. The functions here only restrict the overloads to
// the SoA vector types, so that they do not clash with the overloads of the
// AoS plugins (e.g. Vc AoS) when several plugins are used together.
//
// Note: The plugin has to make the overloads of the math functions (fma,
// sqrt, atan2 and atanh) for its SIMD type visible to argument dependent lookup
// (see algebra::generic::math).
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
  return algebra::generic::math::phi(v);
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
  return algebra::generic::math::perp(v);
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
  return algebra::generic::math::theta(v);
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
/// @returns the cross product of @p a and @p b
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
  requires(N == 3)
DETRAY_HOST_DEVICE constexpr algebra::storage::vector<N, scalar_t, array_t>
cross(const algebra::storage::vector<N, scalar_t, array_t> &a,
      const algebra::storage::vector<N, scalar_t, array_t> &b) {
  return algebra::generic::math::cross(a, b);
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
/// @returns the dot product of @p a and @p b
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr scalar_t dot(
    const algebra::storage::vector<N, scalar_t, array_t> &a,
    const algebra::storage::vector<N, scalar_t, array_t> &b) {
  return algebra::generic::math::dot(a, b);
}

/// Get the norm of a vector in every lane
///
/// @tparam N Number of vector elements
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @param v the input vector
///
/// @returns the euclidean norm of @p v
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr scalar_t norm(
    const algebra::storage::vector<N, scalar_t, array_t> &v) {
  return algebra::generic::math::norm(v);
}

/// Get the normalized vector in every lane
///
/// @tparam N Number of vector elements
/// @tparam scalar_t SIMD scalar type
/// @tparam array_t Underlying storage array template
///
/// @param v the input vector
///
/// @returns the unit vector in the direction of @p v
template <std::size_t N, concepts::simd_scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr algebra::storage::vector<N, scalar_t, array_t>
normalize(const algebra::storage::vector<N, scalar_t, array_t> &v) {
  return algebra::generic::math::normalize(v);
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
  return algebra::generic::math::eta(v);
}

}  // namespace detray::algebra::soa
