// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/algebra/type_traits.hpp"

// System include(s).
#include <concepts>

namespace detray::concepts {

/// Arithmetic types
template <typename T>
concept arithmetic = std::is_arithmetic_v<T>;

template <typename T>
concept arithmetic_cvref = concepts::arithmetic<std::remove_cvref_t<T>>;

/// Enumeration types
template <typename T>
concept enumeration = std::is_enum_v<T>;

// Value concept: Single entry
template <typename T>
concept value = concepts::arithmetic<std::decay_t<T>>;

// Element concept: Single entry, not necessarily for vector/matrix operations
template <typename T>
concept element = concepts::value<T> || concepts::enumeration<std::decay_t<T>>;

/// Scalar concept: Elements of vectors/matrices (can be simd vectors)
template <typename T>
concept scalar =
    !traits::is_matrix<T> && !traits::is_vector<T> && requires(T a, T b) {
      { a + b } -> std::convertible_to<T>;
      { a - b } -> std::convertible_to<T>;
      { a * b } -> std::convertible_to<T>;
      { a / b } -> std::convertible_to<T>;
    };

/// Check if a scalar is simd
template <typename T>
concept simd_scalar = scalar<T> && !std::is_scalar_v<T>;

/// Index concept to access vector/matrix elements
template <typename T>
concept index = std::is_integral_v<T> && !std::same_as<T, bool>;

/// Vector concepts
/// @{
template <typename V>
concept vector = traits::is_vector<V>;

template <typename V>
concept vector2D = vector<V> && (traits::size<V> == 2);

template <typename V>
concept vector3D = vector<V> && (traits::size<V> == 3);
/// @}

/// Point concepts
/// @{
template <typename V>
concept point = vector<V>;

template <typename V>
concept point2D = point<V> && (traits::size<V> == 2);

template <typename V>
concept point3D = point<V> && (traits::size<V> == 3);
/// @}

/// Matrix concepts
/// @{
template <typename M>
concept matrix = traits::is_matrix<M>;

template <typename M>
concept square_matrix = matrix<M> && traits::is_square<M>;

template <typename M>
concept row_matrix = matrix<M> && (traits::rows<M> == 1);

template <typename M>
concept row_matrix3D = row_matrix<M> && (traits::rows<M> == 3);

template <typename M>
concept column_matrix = matrix<M> && (traits::columns<M> == 1);

template <typename M>
concept column_matrix3D = column_matrix<M> && (traits::rows<M> == 3);

template <typename MA, typename MB>
concept matrix_compatible =
    matrix<MA> && matrix<MB> &&
    std::convertible_to<traits::index_t<MA>, traits::index_t<MB>> &&
    std::convertible_to<traits::index_t<MB>, traits::index_t<MA>>;

template <typename MA, typename MB>
concept matrix_multipliable =
    matrix_compatible<MA, MB> && (traits::columns<MA> == traits::rows<MB>) &&
    requires(traits::scalar_t<MA> sa, traits::scalar_t<MB> sb) {
      { (sa * sb) + (sa * sb) };
    };

template <typename MA, typename MB, typename MC>
concept matrix_multipliable_into =
    matrix_multipliable<MA, MB> && matrix_compatible<MA, MC> &&
    matrix_compatible<MB, MC> && (traits::rows<MC> == traits::rows<MA>) &&
    (traits::columns<MC> == traits::columns<MB>) &&
    requires(traits::scalar_t<MA> sa, traits::scalar_t<MB> sb,
             traits::scalar_t<MC>& sc) {
      { sc += (sa * sb) };
    };
/// @}

/// Transform concept
template <typename T>
concept transform3D = requires(T trf) {
  // Local type definitions
  requires scalar<typename T::scalar_type>;
  requires vector3D<typename T::vector3>;
  requires point2D<typename T::point2>;
  requires point3D<typename T::point3>;

  // Methods
  trf.rotation();
  trf.translation();
  trf.point_to_global(typename T::vector3());
  trf.point_to_local(typename T::vector3());
  trf.vector_to_global(typename T::vector3());
  trf.vector_to_local(typename T::vector3());
};

/// Algebra plugin concept
template <typename A>
concept algebra = (concepts::value<typename A::value_type> &&
                   concepts::scalar<typename A::scalar> &&
                   concepts::index<typename A::index_type> &&
                   concepts::vector3D<typename A::vector3D> &&
                   concepts::point2D<typename A::point2D> &&
                   concepts::point3D<typename A::point3D> &&
                   concepts::transform3D<typename A::transform3D> &&
                   concepts::matrix<typename A::template matrix<3, 3>>);

/// Check if an algebra has soa layout
/// @{
template <typename A>
concept soa = (!concepts::arithmetic<get_scalar_t<A>>);

template <typename A>
concept aos = (!concepts::soa<A>);
/// @}

/// Check if a matrix or vector type permits template indexing
/// @{
template <typename M>
concept has_compile_time_2d_access = requires(const M& m) {
  { m.template element<0, 0>() };
};

template <typename V>
concept has_compile_time_1d_access = requires(const V& v) {
  { v.template element<0>() };
};
/// @}
}  // namespace detray::concepts
