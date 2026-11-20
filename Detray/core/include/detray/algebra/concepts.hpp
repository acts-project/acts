// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
///
/// True iff the type is a matrix.
template <typename M>
concept matrix = traits::is_matrix<M>;

/// True iff the type is a matrix of the given shape and scalar type
///
/// Optionally takes a row, column, and scalar
/// type argument; if ROWS is zero, any number of rows is accepted. If COLS is
/// zero, any number of COLS is accepted. If S is void, then any scalar type
/// is accepted. If any of these is a non-default value, the respective
/// quantity or type is statically asserted.
///
/// @note This is defined in terms of @c matrix so that it subsumes it. Every
/// concept below is in turn defined in terms of this one, which keeps the
/// whole matrix hierarchy ordered for overload resolution.
template <typename M, std::size_t ROWS = 0, std::size_t COLS = 0,
          typename S = void>
concept sized_matrix =
    matrix<M> && (ROWS == 0 || traits::rows<M> == ROWS) &&
    (COLS == 0 || traits::columns<M> == COLS) &&
    (std::same_as<S, void> || std::same_as<traits::scalar_t<M>, S>);

template <typename M, std::size_t ROWS = 0, typename S = void>
concept square_matrix = sized_matrix<M, ROWS, ROWS, S> && traits::is_square<M>;

/// True iff the matrix is symmetric
///
/// @note For types of which the symmetry cannot be determined, this is only
/// a readability hint. The concept will be true, even if symmetry can only
/// be determined at runtime.
template <typename M, std::size_t ROWS = 0, typename S = void>
concept symmetric_matrix =
    square_matrix<M, ROWS, S> && (!traits::static_symmetric_assertable_v<M> ||
                                  traits::static_symmetric_v<M>);

/// True iff the matrix is diagonal
///
/// A diagonal matrix is also symmetric, so this refines @c symmetric_matrix.
///
/// @note For types of which the diagonality cannot be determined, this is only
/// a readability hint. The concept will be true, even if diagonality can only
/// be determined at runtime.
template <typename M, std::size_t ROWS = 0, typename S = void>
concept diagonal_matrix =
    symmetric_matrix<M, ROWS, S> &&
    (!traits::static_diagonal_assertable_v<M> || traits::static_diagonal_v<M>);

template <typename M, std::size_t COLS = 0, typename S = void>
concept row_matrix = sized_matrix<M, 1, COLS, S>;

template <typename M, typename S = void>
concept row_matrix3D = row_matrix<M, 3, S>;

template <typename M, std::size_t ROWS = 0, typename S = void>
concept column_matrix = sized_matrix<M, ROWS, 1, S>;

template <typename M, typename S = void>
concept column_matrix3D = column_matrix<M, 3, S>;

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
  requires vector2D<typename T::vector2>;
  requires vector3D<typename T::vector3>;
  requires point2D<typename T::point2>;
  requires point3D<typename T::point3>;

  // Methods
  trf.rotation();
  trf.translation();
  trf.point_to_global(typename T::vector2());
  trf.point_to_global(typename T::vector3());
  trf.point_to_local(typename T::vector3());
  trf.vector_to_global(typename T::vector2());
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
