// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

// Project include(s).
#include "detray/algebra/common/vector.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/type_traits.hpp"

// System include(s).
#include <array>
#include <cassert>

namespace detray::algebra::storage {

/// Generic matrix type that can take vectors as columns
template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t, std::size_t ROW, std::size_t COL>
struct DETRAY_ALIGN(alignof(algebra::storage::vector<ROW, scalar_t, array_t>))
    matrix {
  // The matrix consists of column vectors
  using vector_type = algebra::storage::vector<ROW, scalar_t, array_t>;
  // Value type: Can be simd types
  using scalar_type = scalar_t;

  /// Default constructor
  constexpr matrix() = default;

  /// Construct from given column vectors @param v
  template <concepts::vector... vector_t>
  DETRAY_HOST_DEVICE
    requires(sizeof...(vector_t) == COL)
  explicit matrix(vector_t &&...v) : m_storage{std::forward<vector_t>(v)...} {}

  /// Subscript operator
  /// @{
  DETRAY_HOST_DEVICE
  constexpr const vector_type &operator[](const std::size_t i) const {
    assert(i < COL);
    return m_storage[i];
  }
  DETRAY_HOST_DEVICE
  constexpr vector_type &operator[](const std::size_t i) {
    assert(i < COL);
    return m_storage[i];
  }
  /// @}

  /// @returns the number of rows
  DETRAY_HOST_DEVICE
  static consteval std::size_t rows() { return ROW; }

  /// @returns the number of rows of the underlying vector storage
  /// @note can be different from the matrix rows due to padding
  DETRAY_HOST_DEVICE
  static consteval std::size_t storage_rows() {
    return vector_type::simd_size();
  }

  /// @returns the number of columns
  DETRAY_HOST_DEVICE
  static consteval std::size_t columns() { return COL; }

 private:
  /// Equality operator between two matrices
  template <std::size_t R, std::size_t C, typename S,
            template <typename, std::size_t> class A>
  DETRAY_HOST_DEVICE friend constexpr bool operator==(
      const matrix<A, S, R, C> &lhs, const matrix<A, S, R, C> &rhs);

  /// Sets the trailing uninitialized values to zero.
  /// @{
  // AoS
  template <std::size_t... I>
  DETRAY_HOST_DEVICE
    requires(!std::is_scalar_v<scalar_t>)
  constexpr bool equal(const matrix &rhs, std::index_sequence<I...>) const {
    return (... && (m_storage[I] == rhs[I]));
  }

  // SoA
  template <std::size_t... I>
  DETRAY_HOST
    requires(std::is_scalar_v<scalar_t>)
  constexpr bool equal(const matrix &rhs, std::index_sequence<I...>) const {
    return (... && ((m_storage[I].get() == rhs[I].get()).isFull()));
  }
  /// @}

  /// Arithmetic operators
  /// @{
  template <std::size_t R, std::size_t C, typename S,
            template <typename, std::size_t> class A>
  DETRAY_HOST_DEVICE friend constexpr decltype(auto) operator+(
      const matrix<A, S, R, C> &lhs, const matrix<A, S, R, C> &rhs) noexcept;

  template <std::size_t R, std::size_t C, typename S,
            template <typename, std::size_t> class A>
  DETRAY_HOST_DEVICE friend constexpr decltype(auto) operator-(
      const matrix<A, S, R, C> &lhs, const matrix<A, S, R, C> &rhs) noexcept;

  template <std::size_t R, std::size_t C, typename S1, typename S2,
            template <typename, std::size_t> class A>
  DETRAY_HOST_DEVICE friend constexpr decltype(auto) operator*(
      const S2 a, const matrix<A, S1, R, C> &rhs) noexcept;

  template <std::size_t R, std::size_t C, concepts::scalar S1,
            concepts::scalar S2, template <typename, std::size_t> class A>
  DETRAY_HOST_DEVICE friend constexpr decltype(auto) operator*(
      const matrix<A, S1, R, C> &lhs, const S2 a) noexcept;

  /// Matrix-vector multiplication
  template <std::size_t R, std::size_t C, typename S,
            template <typename, std::size_t> class A>
  DETRAY_HOST_DEVICE friend constexpr decltype(auto) operator*(
      const matrix<A, S, R, C> &lhs, const vector<C, S, A> &v) noexcept;

  /// Matrix-matrix multiplication
  template <std::size_t LR, std::size_t C, std::size_t RC, typename S,
            template <typename, std::size_t> class A>
  DETRAY_HOST_DEVICE friend constexpr decltype(auto) operator*(
      const matrix<A, S, LR, C> &lhs, const matrix<A, S, C, RC> &rhs) noexcept;
  /// @}

  /// Matrix storage
  std::array<vector_type, COL> m_storage;

};  // struct matrix

/// Get a zero-initialized matrix
template <concepts::matrix matrix_t, std::size_t COLS = matrix_t::columns()>
DETRAY_HOST_DEVICE constexpr matrix_t zero() noexcept {
  matrix_t m;

  DETRAY_UNROLL_N(COLS)
  for (std::size_t j = 0u; j < COLS; ++j) {
    // Fill zero initialized vector
    m[j] = typename matrix_t::vector_type{};
  }

  return m;
}

/// Set a matrix to zero
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr void set_zero(matrix_t &m) noexcept {
  m = zero<matrix_t>();
}

/// Build an identity matrix
template <concepts::matrix matrix_t,
          std::size_t R = detray::traits::max_rank<matrix_t>>
DETRAY_HOST_DEVICE constexpr matrix_t identity() noexcept {
  // Zero initialized
  matrix_t m{zero<matrix_t>()};

  DETRAY_UNROLL_N(R)
  for (std::size_t i = 0u; i < R; ++i) {
    m[i][i] = typename matrix_t::scalar_type(1);
  }

  return m;
}

/// Set a matrix to zero
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr void set_identity(matrix_t &m) noexcept {
  m = identity<matrix_t>();
}

/// Transpose the matrix @param m
template <std::size_t ROW, std::size_t COL, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t, std::size_t... I>
DETRAY_HOST_DEVICE constexpr auto transpose(
    const matrix<array_t, scalar_t, ROW, COL> &m,
    std::index_sequence<I...>) noexcept {
  using matrix_T_t = matrix<array_t, scalar_t, COL, ROW>;
  using column_t = typename matrix_T_t::vector_type;

  matrix_T_t res_m;

  DETRAY_UNROLL_N(ROW)
  for (std::size_t j = 0u; j < ROW; ++j) {
    res_m[j] = column_t{m[I][j]...};
  }

  return res_m;
}

/// Build an identity matrix
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr auto transpose(const matrix_t &m) noexcept {
  return transpose(m, std::make_index_sequence<matrix_t::columns()>());
}

/// Equality operator between two matrices
template <std::size_t ROW, std::size_t COL, typename scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr bool operator==(
    const matrix<array_t, scalar_t, ROW, COL> &lhs,
    const matrix<array_t, scalar_t, ROW, COL> &rhs) {
  return lhs.equal(rhs, std::make_index_sequence<COL>());
}

/// Scalar multiplication
template <concepts::matrix matrix_t, concepts::scalar scalar_t,
          std::size_t... J>
DETRAY_HOST_DEVICE constexpr matrix_t matrix_scalar_mul(
    scalar_t a, const matrix_t &rhs, std::index_sequence<J...>) noexcept {
  using mat_scalar_t = detray::traits::scalar_t<matrix_t>;

  return matrix_t{(static_cast<mat_scalar_t>(a) * rhs[J])...};
}

/// Matrix addition
template <concepts::matrix matrix_t, std::size_t... J>
DETRAY_HOST_DEVICE constexpr matrix_t matrix_add(
    const matrix_t &lhs, const matrix_t &rhs,
    std::index_sequence<J...>) noexcept {
  return matrix_t{(lhs[J] + rhs[J])...};
}

template <concepts::matrix matrix_t, std::size_t... J>
DETRAY_HOST_DEVICE constexpr decltype(auto) matrix_sub(
    const matrix_t &lhs, const matrix_t &rhs,
    std::index_sequence<J...>) noexcept {
  return matrix_t{(lhs[J] - rhs[J])...};
}

/// Arithmetic operators
/// @{
template <std::size_t ROW, std::size_t COL, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr decltype(auto) operator+(
    const matrix<array_t, scalar_t, ROW, COL> &lhs,
    const matrix<array_t, scalar_t, ROW, COL> &rhs) noexcept {
  using matrix_t = matrix<array_t, scalar_t, ROW, COL>;

  return matrix_add(lhs, rhs, std::make_index_sequence<matrix_t::columns()>());
}

template <std::size_t ROW, std::size_t COL, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr decltype(auto) operator-(
    const matrix<array_t, scalar_t, ROW, COL> &lhs,
    const matrix<array_t, scalar_t, ROW, COL> &rhs) noexcept {
  using matrix_t = matrix<array_t, scalar_t, ROW, COL>;

  return matrix_sub(lhs, rhs, std::make_index_sequence<matrix_t::columns()>());
}

template <std::size_t R, std::size_t C, typename S1, typename S2,
          template <typename, std::size_t> class A>
DETRAY_HOST_DEVICE constexpr decltype(auto) operator*(
    const S2 a, const matrix<A, S1, R, C> &rhs) noexcept {
  using matrix_t = matrix<A, S2, R, C>;

  return matrix_scalar_mul(static_cast<S1>(a), rhs,
                           std::make_index_sequence<matrix_t::columns()>());
}

template <std::size_t R, std::size_t C, concepts::scalar S1,
          concepts::scalar S2, template <typename, std::size_t> class A>
DETRAY_HOST_DEVICE constexpr decltype(auto) operator*(
    const matrix<A, S1, R, C> &lhs, const S2 a) noexcept {
  return static_cast<S1>(a) * lhs;
}

/// Matrix-vector multiplication
template <std::size_t ROW, std::size_t COL, concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr decltype(auto) operator*(
    const matrix<array_t, scalar_t, ROW, COL> &lhs,
    const vector<COL, scalar_t, array_t> &v) noexcept {
  // Init vector
  vector<ROW, scalar_t, array_t> res_v{v[0] * lhs[0]};

  // Add the rest per column
  DETRAY_UNROLL_N(COL)
  for (std::size_t j = 1u; j < COL; ++j) {
    // fma
    res_v = res_v + v[j] * lhs[j];
  }

  return res_v;
}

/// Matrix-matrix multiplication
template <std::size_t LROW, std::size_t COL, std::size_t RCOL,
          concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr decltype(auto) operator*(
    const matrix<array_t, scalar_t, LROW, COL> &lhs,
    const matrix<array_t, scalar_t, COL, RCOL> &rhs) noexcept {
  matrix<array_t, scalar_t, LROW, RCOL> res_m;

  DETRAY_UNROLL_N(RCOL)
  for (std::size_t j = 0u; j < RCOL; ++j) {
    // Init column j
    res_m[j] = rhs[j][0] * lhs[0];

    // Add the rest per column
    DETRAY_UNROLL_N(COL)
    for (std::size_t i = 1u; i < COL; ++i) {
      // fma
      res_m[j] = res_m[j] + rhs[j][i] * lhs[i];
    }
  }

  return res_m;
}
/// @}

}  // namespace detray::algebra::storage
