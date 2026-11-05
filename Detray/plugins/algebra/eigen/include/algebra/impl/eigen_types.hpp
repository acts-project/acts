// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/detail/eigen_array.hpp"
#include "algebra/impl/eigen_getter.hpp"
#include "algebra/impl/eigen_transform3.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/algebra/type_traits.hpp"

// System include(s).
#include <cstddef>

namespace detray {

namespace algebra::eigen {

/// size type for Eigen storage model
using index_type = int;
/// Value type for Eigen storage model
template <concepts::value T>
using value_type = T;
/// Scalar type for Eigen storage model
template <concepts::value T>
using scalar_type = T;
/// Array type used in the Eigen storage model
template <concepts::scalar T, index_type N>
using storage_type = array<T, N>;
/// Vector type used in the Eigen storage model
template <concepts::scalar T, index_type N>
using vector_type = storage_type<T, N>;
/// Matrix type used in the Eigen storage model
/// If the number of rows is 1, make it RowMajor
template <concepts::scalar T, index_type ROWS, index_type COLS>
using matrix_type =
    Eigen::Matrix<T, ROWS, COLS, static_cast<int>(ROWS == 1), ROWS, COLS>;

/// 3-element "vector" type, using @c eigen::vector_type
template <concepts::scalar T>
using vector3 = storage_type<T, 3>;
/// Point in 3D space, using @c eigen::vector_type
template <concepts::scalar T>
using point3 = vector3<T>;
/// 2-element "vector" type, using @c eigen::vector_type
template <concepts::scalar T>
using vector2 = storage_type<T, 2>;
/// Point in 2D space, using @c eigen::vector_type
template <concepts::scalar T>
using point2 = vector2<T>;

/// Geometry transformation implementation using @c eigen::vector_type
template <concepts::scalar T>
using transform3 = eigen::math::transform3<T>;

/// Element Getter
using element_getter = eigen::storage::element_getter;
/// Block Getter
using block_getter = eigen::storage::block_getter;

}  // namespace algebra::eigen

DETRAY_ALGEBRA_DEFINE_TYPE_TRAITS(algebra::eigen)

// Extra specializations needed for some additional eigen types
namespace traits {

/// Index
/// @{
template <typename derived_t>
struct index<Eigen::MatrixBase<derived_t>> {
  /// Index type used by Eigen algebra traits.
  using type = algebra::eigen::index_type;
};

/// Index trait specialization for a fixed Eigen matrix block.
template <typename T, int ROWS, int COLS, int OPT, int bROWS, int bCOLS>
struct index<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>, bROWS,
                          bCOLS, false>> {
  /// Index type used by Eigen algebra traits.
  using type = algebra::eigen::index_type;
};
/// @}

/// Dimensions
/// @{
/// Dimensions trait for a fixed Eigen matrix expression.
template <typename derived_t>
struct dimensions<Eigen::MatrixBase<derived_t>> {
  /// Index type used for compile-time dimensions.
  using index_type = index_t<Eigen::MatrixBase<derived_t>>;

  /// Number of matrix dimensions.
  static constexpr index_type _dim{2};
  /// Compile-time row count.
  static constexpr index_type _rows{
      Eigen::MatrixBase<derived_t>::RowsAtCompileTime};
  /// Compile-time column count.
  static constexpr index_type _columns{
      Eigen::MatrixBase<derived_t>::ColsAtCompileTime};
};

/// Dimensions trait for a fixed-size Eigen matrix block.
template <typename T, int ROWS, int COLS, int OPT, int bROWS, int bCOLS>
struct dimensions<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>,
                               bROWS, bCOLS, false>> {
  /// Index type used for compile-time dimensions.
  using index_type =
      index_t<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>, bROWS,
                           bCOLS, false>>;

  /// Number of matrix dimensions.
  static constexpr index_type _dim{2};
  /// Block row count.
  static constexpr index_type _rows{bROWS};
  /// Block column count.
  static constexpr index_type _columns{bCOLS};
};
/// @}

/// Value
/// @{
template <typename derived_t>
struct value<Eigen::MatrixBase<derived_t>> {
  /// Scalar element type of the Eigen expression.
  using type = typename Eigen::MatrixBase<derived_t>::value_type;
};

/// Value trait specialization for a fixed Eigen matrix block.
template <typename T, int ROWS, int COLS, int OPT, int bROWS, int bCOLS>
struct value<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>, bROWS,
                          bCOLS, false>> {
  /// Scalar element type stored by the block.
  using type = T;
};
/// @}

/// Vector
/// @{
/// Vector trait for an Eigen expression with a single compile-time column.
template <typename derived_t>
struct vector<Eigen::MatrixBase<derived_t>> {
  /// Construct a vector representation with another scalar and length.
  template <typename other_T, int other_N>
  using other_type = algebra::eigen::vector_type<other_T, other_N>;

  /// Vector type matching the expression's scalar type and row count.
  using type = other_type<value_t<Eigen::MatrixBase<derived_t>>,
                          rows<Eigen::MatrixBase<derived_t>>>;
};
/// @}

/// Matrix
/// @{
/// Matrix trait for an Eigen matrix expression.
template <typename derived_t>
struct matrix<Eigen::MatrixBase<derived_t>> {
  /// Construct a matrix representation with the requested type and dimensions.
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type = algebra::eigen::matrix_type<
      typename Eigen::MatrixBase<derived_t>::value_type,
      Eigen::MatrixBase<derived_t>::RowsAtCompileTime,
      Eigen::MatrixBase<derived_t>::ColsAtCompileTime>;

  /// Eigen matrix expression type itself.
  using type = Eigen::MatrixBase<derived_t>;
};

/// Matrix trait for an Eigen block expression.
template <typename T, int ROWS, int COLS, int OPT, int bROWS, int bCOLS>
struct matrix<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>, bROWS,
                           bCOLS, false>> {
  /// Construct an Eigen matrix type with the requested dimensions.
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type =
      algebra::eigen::matrix_type<other_T, other_ROWS, other_COLS>;

  /// Eigen block type represented by this trait.
  using type = Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>,
                            bROWS, bCOLS, false>;
};
/// @}

/// Element getter
/// @{
template <typename derived_t>
struct element_getter<Eigen::MatrixBase<derived_t>> {
  /// Accessor type for Eigen matrix elements.
  using type = algebra::eigen::storage::element_getter;
};

/// Element accessor trait for an Eigen matrix block.
template <typename T, int ROWS, int COLS, int OPT, int bROWS, int bCOLS>
struct element_getter<Eigen::Block<
    Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>, bROWS, bCOLS, false>> {
  /// Accessor type for Eigen matrix elements.
  using type = algebra::eigen::storage::element_getter;
};
/// @}

/// Block getter
/// @{
template <typename derived_t>
struct block_getter<Eigen::MatrixBase<derived_t>> {
  /// Accessor type for Eigen matrix blocks.
  using type = algebra::eigen::storage::block_getter;
};
/// @}

}  // namespace traits

}  // namespace detray
