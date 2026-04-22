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
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/assert.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s).
#include <cstddef>
#include <type_traits>

namespace detray::algebra::array::storage {

/// "Element getter", assuming a simple 2D array access
struct element_getter {
  /// Operator getting a reference to one element of a non-const matrix
  template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  DETRAY_HOST_DEVICE constexpr scalar_t &operator()(
      array_t<array_t<scalar_t, ROWS>, COLS> &m, std::size_t row,
      std::size_t col) const {
    assert(row < ROWS);
    assert(col < COLS);
    return m[col][row];
  }

  /// Operator getting one value of a const matrix
  template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  DETRAY_HOST_DEVICE constexpr scalar_t operator()(
      const array_t<array_t<scalar_t, ROWS>, COLS> &m, std::size_t row,
      std::size_t col) const {
    assert(row < ROWS);
    assert(col < COLS);
    return m[col][row];
  }

  /// Operator getting a reference to one element of a non-const matrix
  template <std::size_t N, concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  DETRAY_HOST_DEVICE constexpr scalar_t &operator()(
      array_t<array_t<scalar_t, N>, 1> &m, std::size_t row) const {
    assert(row < N);
    return m[0][row];
  }

  /// Operator getting a reference to one element of a const matrix
  template <std::size_t N, concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  DETRAY_HOST_DEVICE constexpr scalar_t operator()(
      const array_t<array_t<scalar_t, N>, 1> &m, std::size_t row) const {
    assert(row < N);
    return m[0][row];
  }

  /// Operator getting a reference to one element of a non-const matrix
  template <std::size_t N, concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  DETRAY_HOST_DEVICE constexpr scalar_t &operator()(array_t<scalar_t, N> &m,
                                                    std::size_t row) const {
    assert(row < N);
    return m[row];
  }

  /// Operator getting a reference to one element of a const matrix
  template <std::size_t N, concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  DETRAY_HOST_DEVICE constexpr scalar_t operator()(
      const array_t<scalar_t, N> &m, std::size_t row) const {
    assert(row < N);
    return m[row];
  }

};  // struct element_getter

/// Function extracting an element from a matrix (const)
template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t, std::size_t ROWS, std::size_t COLS>
DETRAY_HOST_DEVICE constexpr scalar_t element(
    const array_t<array_t<scalar_t, ROWS>, COLS> &m, std::size_t row,
    std::size_t col) {
  return element_getter()(m, row, col);
}

/// Function extracting an element from a matrix (non-const)
template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t, std::size_t ROWS, std::size_t COLS>
DETRAY_HOST_DEVICE constexpr scalar_t &element(
    array_t<array_t<scalar_t, ROWS>, COLS> &m, std::size_t row,
    std::size_t col) {
  return element_getter()(m, row, col);
}

/// Function extracting an element from a matrix (const)
template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t, std::size_t N>
DETRAY_HOST_DEVICE constexpr scalar_t element(
    const array_t<array_t<scalar_t, N>, 1> &m, std::size_t row) {
  return element_getter()(m, row);
}

/// Function extracting an element from a matrix (non-const)
template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t, std::size_t N>
DETRAY_HOST_DEVICE constexpr scalar_t &element(
    array_t<array_t<scalar_t, N>, 1> &m, std::size_t row) {
  return element_getter()(m, row);
}

/// Function extracting an element from a vector (const)
template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t, std::size_t N>
DETRAY_HOST_DEVICE constexpr scalar_t element(const array_t<scalar_t, N> &v,
                                              std::size_t row) {
  return element_getter()(v, row);
}

/// Function extracting an element from a vector (non-const)
template <template <typename, std::size_t> class array_t,
          concepts::scalar scalar_t, std::size_t N>
DETRAY_HOST_DEVICE constexpr scalar_t &element(array_t<scalar_t, N> &v,
                                               std::size_t row) {
  return element_getter()(v, row);
}

template <std::size_t I, std::size_t J, concepts::matrix M>
DETRAY_HOST_DEVICE decltype(auto) element(M &matrix) {
  if constexpr (concepts::has_compile_time_2d_access<M>) {
    return matrix.template element<I, J>();
  } else {
    using index_t = detray::traits::index_t<std::decay_t<M>>;
    return element(matrix, index_t(I), index_t(J));
  }
}

template <std::size_t I, concepts::vector V>
DETRAY_HOST_DEVICE decltype(auto) element(V &vector) {
  if constexpr (concepts::has_compile_time_1d_access<V>) {
    return vector.template element<I>();
  } else {
    using index_t = detray::traits::index_t<std::decay_t<V>>;
    return element(vector, index_t(I));
  }
}

/// "Block getter", assuming a simple 2D array access
struct block_getter {
  /// Operator producing a sub-matrix from a const matrix
  template <std::size_t ROWS, std::size_t COLS, std::size_t oROWS,
            std::size_t oCOLS, concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  DETRAY_HOST_DEVICE constexpr auto operator()(
      const array_t<array_t<scalar_t, oROWS>, oCOLS> &m, std::size_t row,
      std::size_t col) const {
    array_t<array_t<scalar_t, ROWS>, COLS> submatrix{};

    for (std::size_t icol = col; icol < col + COLS; ++icol) {
      for (std::size_t irow = row; irow < row + ROWS; ++irow) {
        submatrix[icol - col][irow - row] = m[icol][irow];
      }
    }
    return submatrix;
  }

  /// Operator producing a vector out of a const matrix
  template <std::size_t SIZE, std::size_t ROWS, std::size_t COLS,
            concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  DETRAY_HOST_DEVICE constexpr array_t<scalar_t, SIZE> vector(
      const array_t<array_t<scalar_t, ROWS>, COLS> &m, std::size_t row,
      std::size_t col) const {
    assert(col < COLS);
    assert(row + SIZE <= ROWS);

    array_t<scalar_t, SIZE> subvector{};

    for (std::size_t irow = row; irow < row + SIZE; ++irow) {
      subvector[irow - row] = m[col][irow];
    }

    return subvector;
  }

  /// Sets a matrix of dimension @tparam ROW and @tparam COL as submatrix of
  /// @param m beginning at row @param row and column @param col
  template <std::size_t ROWS, std::size_t COLS, class input_matrix_type,
            concepts::scalar scalar_t,
            template <typename, std::size_t> class array_t>
  DETRAY_HOST_DEVICE constexpr void set(
      input_matrix_type &m, const array_t<array_t<scalar_t, ROWS>, COLS> &b,
      std::size_t row, std::size_t col) const {
    for (std::size_t j = 0u; j < COLS; ++j) {
      for (std::size_t i = 0u; i < ROWS; ++i) {
        m[j + col][i + row] = b[j][i];
      }
    }
  }

  /// Sets a vector of length @tparam ROW as submatrix of
  /// @param m beginning at row @param row and column @param col
  template <std::size_t ROWS, concepts::scalar scalar_t,
            template <typename, std::size_t> class vector_t,
            class input_matrix_type>
  DETRAY_HOST_DEVICE constexpr void set(input_matrix_type &m,
                                        const vector_t<scalar_t, ROWS> &b,
                                        std::size_t row,
                                        std::size_t col) const {
    for (std::size_t i = 0; i < ROWS; ++i) {
      m[col][i + row] = b[i];
    }
  }

};  // struct block_getter

/// @returns a matrix of dimension @tparam ROW and @tparam COL that contains
/// the submatrix of @param m beginning at row @param row and column
/// @param col
template <std::size_t ROWS, std::size_t COLS, class input_matrix_type>
DETRAY_HOST_DEVICE constexpr decltype(auto) block(const input_matrix_type &m,
                                                  std::size_t row,
                                                  std::size_t col) {
  constexpr block_getter block{};
  return block.template operator()<ROWS, COLS>(m, row, col);
}

/// Function extracting a vector from a matrix
template <std::size_t SIZE, std::size_t ROWS, std::size_t COLS,
          concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr array_t<scalar_t, SIZE> vector(
    const array_t<array_t<scalar_t, ROWS>, COLS> &m, std::size_t row,
    std::size_t col) {
  constexpr block_getter block{};
  return block.template vector<SIZE>(m, row, col);
}

/// Sets a matrix of dimension @tparam ROW and @tparam COL as submatrix of
/// @param m beginning at row @param row and column @param col
template <std::size_t ROWS, std::size_t COLS, class input_matrix_type,
          concepts::scalar scalar_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE constexpr void set_block(
    input_matrix_type &m, const array_t<array_t<scalar_t, ROWS>, COLS> &b,
    std::size_t row, std::size_t col) {
  constexpr block_getter block{};
  block.set(m, b, row, col);
}

/// Sets a vector of length @tparam ROW as submatrix of
/// @param m beginning at row @param row and column @param col
template <std::size_t ROWS, concepts::scalar scalar_t,
          template <typename, std::size_t> class vector_t,
          class input_matrix_type>
DETRAY_HOST_DEVICE constexpr void set_block(input_matrix_type &m,
                                            const vector_t<scalar_t, ROWS> &b,
                                            std::size_t row, std::size_t col) {
  constexpr block_getter block{};
  block.set(m, b, row, col);
}

}  // namespace detray::algebra::array::storage
