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
#include "detray/definitions/detail/qualifiers.hpp"

// Eigen include(s).
#ifdef _MSC_VER
#pragma warning(push, 0)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 20012
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__
#include <Eigen/Core>
#ifdef _MSC_VER
#pragma warning(pop)
#endif  // MSVC
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#endif  // __NVCC_DIAG_PRAGMA_SUPPORT__

// System include(s).
#include <type_traits>

namespace detray::algebra::eigen::storage {

/// Functor used to access elements of Eigen matrices
struct element_getter {
  /// Get non-const access to a matrix element
  template <typename derived_type, concepts::index index_t_1,
            concepts::index index_t_2>
    requires std::is_base_of_v<
        Eigen::DenseCoeffsBase<derived_type, Eigen::WriteAccessors>,
        Eigen::MatrixBase<derived_type>>
  DETRAY_HOST_DEVICE constexpr auto &operator()(
      Eigen::MatrixBase<derived_type> &m, index_t_1 row, index_t_2 col) const {
    return m(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
  }
  /// Get const access to a matrix element
  template <typename derived_type, concepts::index index_t_1,
            concepts::index index_t_2>
  DETRAY_HOST_DEVICE constexpr auto operator()(
      const Eigen::MatrixBase<derived_type> &m, index_t_1 row,
      index_t_2 col) const {
    return m(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
  }
  /// Get non-const access to a matrix element
  template <typename derived_type, concepts::index index_t>
    requires std::is_base_of_v<
        Eigen::DenseCoeffsBase<derived_type, Eigen::WriteAccessors>,
        Eigen::MatrixBase<derived_type>>
  DETRAY_HOST_DEVICE constexpr auto &operator()(
      Eigen::MatrixBase<derived_type> &m, index_t row) const {
    return m(static_cast<Eigen::Index>(row));
  }
  /// Get const access to a matrix element
  template <typename derived_type, concepts::index index_t>
  DETRAY_HOST_DEVICE constexpr auto operator()(
      const Eigen::MatrixBase<derived_type> &m, index_t row) const {
    return m(static_cast<Eigen::Index>(row));
  }
};  // struct element_getter

/// Function extracting an element from a matrix (const)
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr decltype(auto) element(
    const Eigen::MatrixBase<derived_type> &m, std::size_t row,
    std::size_t col) {
  return element_getter()(m, static_cast<Eigen::Index>(row),
                          static_cast<Eigen::Index>(col));
}

/// Function extracting an element from a matrix (non-const)
template <typename derived_type>
  requires std::is_base_of_v<
      Eigen::DenseCoeffsBase<derived_type, Eigen::WriteAccessors>,
      Eigen::MatrixBase<derived_type>>
DETRAY_HOST_DEVICE constexpr decltype(auto) element(
    Eigen::MatrixBase<derived_type> &m, std::size_t row, std::size_t col) {
  return element_getter()(m, static_cast<Eigen::Index>(row),
                          static_cast<Eigen::Index>(col));
}

/// Function extracting an element from a matrix (const)
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr decltype(auto) element(
    const Eigen::MatrixBase<derived_type> &m, std::size_t row) {
  return element_getter()(m, static_cast<Eigen::Index>(row));
}

/// Function extracting an element from a matrix (non-const)
template <typename derived_type>
  requires std::is_base_of_v<
      Eigen::DenseCoeffsBase<derived_type, Eigen::WriteAccessors>,
      Eigen::MatrixBase<derived_type>>
DETRAY_HOST_DEVICE constexpr decltype(auto) element(
    Eigen::MatrixBase<derived_type> &m, std::size_t row) {
  return element_getter()(m, static_cast<Eigen::Index>(row));
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

/// Functor used to extract a block from Eigen matrices
struct block_getter {
  template <int kROWS, int kCOLS, typename derived_type,
            concepts::index index_t_1, concepts::index index_t_2>
  DETRAY_HOST_DEVICE decltype(auto) operator()(
      const Eigen::MatrixBase<derived_type> &m, const index_t_1 row,
      const index_t_2 col) const {
    return m.template block<kROWS, kCOLS>(row, col);
  }

  template <int kROWS, int kCOLS, typename derived_type,
            concepts::index index_t_1, concepts::index index_t_2>
  DETRAY_HOST_DEVICE decltype(auto) operator()(
      Eigen::MatrixBase<derived_type> &m, const index_t_1 row,
      const index_t_2 col) const {
    return m.template block<kROWS, kCOLS>(row, col);
  }

  template <int SIZE, typename derived_type, concepts::index index_t_1,
            concepts::index index_t_2>
  DETRAY_HOST_DEVICE decltype(auto) vector(Eigen::MatrixBase<derived_type> &m,
                                           const index_t_1 row,
                                           const index_t_2 col) const {
    return m.template block<SIZE, 1>(row, col);
  }

  template <int SIZE, typename derived_type, concepts::index index_t_1,
            concepts::index index_t_2>
  DETRAY_HOST_DEVICE decltype(auto) vector(
      const Eigen::MatrixBase<derived_type> &m, const index_t_1 row,
      const index_t_2 col) const {
    return m.template block<SIZE, 1>(row, col);
  }

  /// Operator setting a block
  template <typename derived_type1, typename derived_type2,
            concepts::index index_t_1, concepts::index index_t_2>
  DETRAY_HOST_DEVICE void set(Eigen::MatrixBase<derived_type1> &m,
                              const Eigen::MatrixBase<derived_type2> &b,
                              const index_t_1 row, const index_t_2 col) const {
    using block_t = Eigen::MatrixBase<derived_type2>;
    constexpr auto R{block_t::RowsAtCompileTime};
    constexpr auto C{block_t::ColsAtCompileTime};
    m.template block<R, C>(static_cast<Eigen::Index>(row),
                           static_cast<Eigen::Index>(col)) = b;
  }

};  // struct block_getter

/// Operator getting a block of a const matrix
template <int ROWS, int COLS, class derived_type>
DETRAY_HOST_DEVICE decltype(auto) block(
    const Eigen::MatrixBase<derived_type> &m, std::size_t row,
    std::size_t col) {
  return block_getter{}.template operator()<ROWS, COLS>(
      m, static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
}

/// Operator getting a block of a non-const matrix
template <int ROWS, int COLS, class derived_type>
DETRAY_HOST_DEVICE decltype(auto) block(Eigen::MatrixBase<derived_type> &m,
                                        std::size_t row, std::size_t col) {
  return block_getter{}.template operator()<ROWS, COLS>(
      m, static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
}

/// Function extracting a slice from the matrix
template <int SIZE, typename derived_type>
DETRAY_HOST_DEVICE constexpr decltype(auto) vector(
    const Eigen::MatrixBase<derived_type> &m, std::size_t row,
    std::size_t col) {
  return block_getter{}.template vector<SIZE>(m, static_cast<Eigen::Index>(row),
                                              static_cast<Eigen::Index>(col));
}

/// Operator setting a block
template <typename derived_type1, typename derived_type2>
DETRAY_HOST_DEVICE void set_block(Eigen::MatrixBase<derived_type1> &m,
                                  const Eigen::MatrixBase<derived_type2> &b,
                                  std::size_t row, std::size_t col) {
  block_getter{}.set(m, b, static_cast<Eigen::Index>(row),
                     static_cast<Eigen::Index>(col));
}

}  // namespace detray::algebra::eigen::storage
