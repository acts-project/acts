// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// ROOT/Smatrix include(s).
#include <Math/Expression.h>
#include <Math/Functions.h>
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TMath.h>

// System include(s).
#include <cassert>

namespace detray::algebra::smatrix::storage {

/// Functor used to access elements of SMatrix matrices
struct element_getter {
  template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t &operator()(
      ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, unsigned int row,
      unsigned int col) const {
    assert(row < ROWS);
    assert(col < COLS);
    return m(row, col);
  }

  template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t operator()(
      const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, unsigned int row,
      unsigned int col) const {
    assert(row < ROWS);
    assert(col < COLS);
    return m(row, col);
  }

  template <unsigned int N, concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t &operator()(
      ROOT::Math::SMatrix<scalar_t, N, 1> &m, unsigned int row) const {
    assert(row < N);
    return m(row);
  }

  template <unsigned int N, concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t operator()(
      const ROOT::Math::SMatrix<scalar_t, N, 1> &m, unsigned int row) const {
    assert(row < N);
    return m(row, 0);
  }

  template <concepts::scalar scalar_t, unsigned int N>
  DETRAY_HOST_DEVICE constexpr scalar_t &operator()(
      ROOT::Math::SVector<scalar_t, N> &m, unsigned int row) const {
    assert(row < N);
    return m(row);
  }

  template <unsigned int N, concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t operator()(
      const ROOT::Math::SVector<scalar_t, N> &m, unsigned int row) const {
    assert(row < N);
    return m(row);
  }
};  // element_getter

/// Function extracting an element from a matrix (const)
template <concepts::scalar scalar_t, unsigned int ROWS, unsigned int COLS>
DETRAY_HOST_DEVICE constexpr scalar_t element(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, std::size_t row,
    std::size_t col) {
  return element_getter()(m, static_cast<unsigned int>(row),
                          static_cast<unsigned int>(col));
}

/// Function extracting an element from a matrix (non-const)
template <concepts::scalar scalar_t, unsigned int ROWS, unsigned int COLS>
DETRAY_HOST_DEVICE constexpr scalar_t &element(
    ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, std::size_t row,
    std::size_t col) {
  return element_getter()(m, static_cast<unsigned int>(row),
                          static_cast<unsigned int>(col));
}

/// Function extracting an element from a matrix (const)
template <concepts::scalar scalar_t, unsigned int N>
DETRAY_HOST_DEVICE constexpr scalar_t element(
    const ROOT::Math::SMatrix<scalar_t, N, 1> &m, std::size_t row) {
  return element_getter()(m, static_cast<unsigned int>(row));
}

/// Function extracting an element from a matrix (non-const)
template <concepts::scalar scalar_t, unsigned int N>
DETRAY_HOST_DEVICE constexpr scalar_t &element(
    ROOT::Math::SMatrix<scalar_t, N, 1> &m, std::size_t row) {
  return element_getter()(m, static_cast<unsigned int>(row));
}

/// Function extracting an element from a matrix (const)
template <concepts::scalar scalar_t, unsigned int N>
DETRAY_HOST_DEVICE constexpr scalar_t element(
    const ROOT::Math::SVector<scalar_t, N> &m, std::size_t row) {
  return element_getter()(m, static_cast<unsigned int>(row));
}

/// Function extracting an element from a matrix (non-const)
template <concepts::scalar scalar_t, unsigned int N>
DETRAY_HOST_DEVICE constexpr scalar_t &element(
    ROOT::Math::SVector<scalar_t, N> &m, std::size_t row) {
  return element_getter()(m, static_cast<unsigned int>(row));
}

template <std::size_t I, std::size_t J, concepts::matrix M>
DETRAY_HOST_DEVICE decltype(auto) element(M &matrix) {
  if constexpr (concepts::has_compile_time_2d_access<M>) {
    return matrix.template element<I, J>();
  } else {
    using index_t = detray::traits::index_t<std::decay_t<M>>;
    return element(matrix, static_cast<index_t>(I), static_cast<index_t>(J));
  }
}

template <std::size_t I, concepts::vector V>
DETRAY_HOST_DEVICE decltype(auto) element(V &vector) {
  if constexpr (concepts::has_compile_time_1d_access<V>) {
    return vector.template element<I>();
  } else {
    using index_t = detray::traits::index_t<std::decay_t<V>>;
    return element(vector, static_cast<index_t>(I));
  }
}

/// Functor used to extract a block from SMatrix matrices
struct block_getter {
  template <unsigned int ROWS, unsigned int COLS, unsigned int oROWS,
            unsigned int oCOLS, concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE ROOT::Math::SMatrix<scalar_t, ROWS, COLS> operator()(
      const ROOT::Math::SMatrix<scalar_t, oROWS, oCOLS> &m, unsigned int row,
      unsigned int col) const {
    return m.template Sub<ROOT::Math::SMatrix<scalar_t, ROWS, COLS>>(row, col);
  }

  template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE ROOT::Math::SVector<scalar_t, SIZE> vector(
      const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, unsigned int row,
      unsigned int col) const {
    // TODO: SMatrix bug?
    // return m.template SubCol<ROOT::Math::SVector<scalar_t, SIZE>>(col,
    // row);

    assert(col < COLS);
    assert(row + SIZE <= ROWS);

    ROOT::Math::SVector<scalar_t, SIZE> ret;

    for (unsigned int irow = row; irow < row + SIZE; ++irow) {
      ret[irow - row] = m(irow, col);
    }

    return ret;
  }

  /// Operator setting a block with a matrix
  template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t,
            concepts::matrix input_matrix_type>
  DETRAY_HOST_DEVICE void set(
      input_matrix_type &m, const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &b,
      std::size_t row, std::size_t col) const {
    for (unsigned int i = 0; i < ROWS; ++i) {
      for (unsigned int j = 0; j < COLS; ++j) {
        m(i + static_cast<unsigned int>(row),
          j + static_cast<unsigned int>(col)) = b(i, j);
      }
    }
  }

  /// Operator setting a block with a vector
  template <unsigned int ROWS, concepts::scalar scalar_t,
            concepts::matrix input_matrix_type>
  DETRAY_HOST_DEVICE void set(input_matrix_type &m,
                              const ROOT::Math::SVector<scalar_t, ROWS> &b,
                              unsigned int row, unsigned int col) const {
    for (unsigned int i = 0; i < ROWS; ++i) {
      m(i + row, col) = b[i];
    }
  }

};  // struct block_getter

/// Operator getting a block of a const matrix
template <unsigned int ROWS, unsigned int COLS, unsigned int oROWS,
          unsigned int oCOLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE ROOT::Math::SMatrix<scalar_t, ROWS, COLS> block(
    const ROOT::Math::SMatrix<scalar_t, oROWS, oCOLS> &m, std::size_t row,
    std::size_t col) {
  return block_getter{}.template operator()<ROWS, COLS>(
      m, static_cast<unsigned int>(row), static_cast<unsigned int>(col));
}

/// Function extracting a slice from the matrix used by
/// @c algebra::smatrix::transform3
template <unsigned int SIZE, unsigned int ROWS, unsigned int COLS,
          concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr auto vector(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m, std::size_t row,
    std::size_t col) {
  return block_getter{}.template vector<SIZE>(m, static_cast<unsigned int>(row),
                                              static_cast<unsigned int>(col));
}

/// Operator setting a block with a matrix
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t,
          concepts::matrix input_matrix_type>
DETRAY_HOST_DEVICE void set_block(
    input_matrix_type &m, const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &b,
    std::size_t row, std::size_t col) {
  block_getter{}.set(m, b, row, col);
}

/// Operator setting a block with a vector
template <unsigned int ROWS, concepts::scalar scalar_t,
          concepts::matrix input_matrix_type>
DETRAY_HOST_DEVICE void set_block(input_matrix_type &m,
                                  const ROOT::Math::SVector<scalar_t, ROWS> &b,
                                  unsigned int row, unsigned int col) {
  block_getter{}.set(m, b, row, col);
}

}  // namespace detray::algebra::smatrix::storage
