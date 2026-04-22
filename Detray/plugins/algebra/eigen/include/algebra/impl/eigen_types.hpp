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
  using type = algebra::eigen::index_type;
};

template <typename T, int ROWS, int COLS, int OPT, int bROWS, int bCOLS>
struct index<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>, bROWS,
                          bCOLS, false>> {
  using type = algebra::eigen::index_type;
};
/// @}

/// Dimensions
/// @{
template <typename derived_t>
struct dimensions<Eigen::MatrixBase<derived_t>> {
  using index_type = index_t<Eigen::MatrixBase<derived_t>>;

  static constexpr index_type _dim{2};
  static constexpr index_type _rows{
      Eigen::MatrixBase<derived_t>::RowsAtCompileTime};
  static constexpr index_type _columns{
      Eigen::MatrixBase<derived_t>::ColsAtCompileTime};
};

template <typename T, int ROWS, int COLS, int OPT, int bROWS, int bCOLS>
struct dimensions<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>,
                               bROWS, bCOLS, false>> {
  using index_type =
      index_t<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>, bROWS,
                           bCOLS, false>>;

  static constexpr index_type _dim{2};
  static constexpr index_type _rows{bROWS};
  static constexpr index_type _columns{bCOLS};
};
/// @}

/// Value
/// @{
template <typename derived_t>
struct value<Eigen::MatrixBase<derived_t>> {
  using type = typename Eigen::MatrixBase<derived_t>::value_type;
};

template <typename T, int ROWS, int COLS, int OPT, int bROWS, int bCOLS>
struct value<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>, bROWS,
                          bCOLS, false>> {
  using type = T;
};
/// @}

/// Vector
/// @{
template <typename derived_t>
struct vector<Eigen::MatrixBase<derived_t>> {
  template <typename other_T, int other_N>
  using other_type = algebra::eigen::vector_type<other_T, other_N>;

  using type = other_type<value_t<Eigen::MatrixBase<derived_t>>,
                          rows<Eigen::MatrixBase<derived_t>>>;
};
/// @}

/// Matrix
/// @{
template <typename derived_t>
struct matrix<Eigen::MatrixBase<derived_t>> {
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type = algebra::eigen::matrix_type<
      typename Eigen::MatrixBase<derived_t>::value_type,
      Eigen::MatrixBase<derived_t>::RowsAtCompileTime,
      Eigen::MatrixBase<derived_t>::ColsAtCompileTime>;

  using type = Eigen::MatrixBase<derived_t>;
};

template <typename T, int ROWS, int COLS, int OPT, int bROWS, int bCOLS>
struct matrix<Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>, bROWS,
                           bCOLS, false>> {
  template <typename other_T, int other_ROWS, int other_COLS>
  using other_type =
      algebra::eigen::matrix_type<other_T, other_ROWS, other_COLS>;

  using type = Eigen::Block<Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>,
                            bROWS, bCOLS, false>;
};
/// @}

/// Element getter
/// @{
template <typename derived_t>
struct element_getter<Eigen::MatrixBase<derived_t>> {
  using type = algebra::eigen::storage::element_getter;
};

template <typename T, int ROWS, int COLS, int OPT, int bROWS, int bCOLS>
struct element_getter<Eigen::Block<
    Eigen::Matrix<T, ROWS, COLS, OPT, ROWS, COLS>, bROWS, bCOLS, false>> {
  using type = algebra::eigen::storage::element_getter;
};
/// @}

/// Block getter
/// @{
template <typename derived_t>
struct block_getter<Eigen::MatrixBase<derived_t>> {
  using type = algebra::eigen::storage::block_getter;
};
/// @}

}  // namespace traits

}  // namespace detray
