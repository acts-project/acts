// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "algebra/impl/fastor_types.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// Fastor include(s).
#ifdef _MSC_VER
#pragma warning(disable : 4244 4701 4702)
#endif  // MSVC
#include <Fastor/Fastor.h>
#ifdef _MSC_VER
#pragma warning(default : 4244 4701 4702)
#endif  // MSVC

namespace detray::algebra::fastor::math {

/// Create zero matrix
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr matrix_t zero() {
  return matrix_t(0);
}

/// Create identity matrix
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr matrix_t identity() {
  using scalar_t = detray::traits::value_t<matrix_t>;
  constexpr auto rows{detray::traits::rows<matrix_t>};
  constexpr auto cols{detray::traits::columns<matrix_t>};

  if constexpr (rows == cols) {
    matrix_type<scalar_t, rows, cols> identity_matrix;
    identity_matrix.eye2();
    return identity_matrix;
  } else if constexpr (rows > cols) {
    matrix_type<scalar_t, rows, rows> identity_matrix;
    identity_matrix.eye2();
    return matrix_t(
        identity_matrix(Fastor::fseq<0, rows>(), Fastor::fseq<0, cols>()));
  } else {
    matrix_type<scalar_t, cols, cols> identity_matrix;
    identity_matrix.eye2();
    return matrix_t(
        identity_matrix(Fastor::fseq<0, rows>(), Fastor::fseq<0, cols>()));
  }
}

/// Set input matrix as zero matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr void set_zero(
    matrix_type<scalar_t, ROWS, COLS> &m) {
  m.zeros();
}

/// Set input matrix as identity matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr void set_identity(
    matrix_type<scalar_t, ROWS, COLS> &m) {
  m = identity<matrix_type<scalar_t, ROWS, COLS>>();
}

/// Create transpose matrix
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS> transpose(
    const matrix_type<scalar_t, ROWS, COLS> &m) {
  return Fastor::transpose(m);
}

/// Column-wise cross product
/// @{
template <typename derived_1_t, typename derived_2_t>
DETRAY_HOST_DEVICE constexpr auto column_wise_cross(
    const Fastor::AbstractTensor<derived_1_t, 2> &m,
    const Fastor::AbstractTensor<derived_2_t, 1> &v) {
  using scalar_t = typename derived_1_t::scalar_type;

  matrix_type<scalar_t, 3, 3> ret;
  for (std::size_t j = 0u; j < 3; j++) {
    ret(Fastor::all, j) =
        Fastor::cross(static_cast<vector_type<scalar_t, 3>>(
                          Fastor::evaluate(m)(Fastor::all, j)),
                      Fastor::evaluate(v));
  }

  return ret;
}

template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr auto column_wise_cross(
    const matrix_type<scalar_t, COLS, ROWS> &m,
    const vector_type<scalar_t, ROWS> &v) {
  matrix_type<scalar_t, ROWS, ROWS> ret;
  for (std::size_t j = 0u; j < ROWS; j++) {
    ret(Fastor::all, j) = Fastor::cross(
        static_cast<vector_type<scalar_t, ROWS>>(m(Fastor::all, j)), v);
  }

  return ret;
}
/// @}

/// Column-wise product
/// @{
template <typename derived_1_t, typename derived_2_t>
DETRAY_HOST_DEVICE constexpr auto column_wise_multiply(
    const Fastor::AbstractTensor<derived_1_t, 2> &m,
    const Fastor::AbstractTensor<derived_2_t, 1> &v) {
  using scalar_t = typename derived_1_t::scalar_type;
  constexpr std::size_t N{derived_1_t::size()};

  matrix_type<scalar_t, N, N> ret;
  for (std::size_t j = 0u; j < N; j++) {
    ret(Fastor::all, j) = m.self()(Fastor::all, j) * v;
  }

  return ret;
}

template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr auto column_wise_multiply(
    const matrix_type<scalar_t, COLS, ROWS> &m,
    const vector_type<scalar_t, ROWS> &v) {
  matrix_type<scalar_t, ROWS, ROWS> ret;
  for (std::size_t j = 0u; j < ROWS; j++) {
    ret(Fastor::all, j) =
        static_cast<vector_type<scalar_t, ROWS>>(m(Fastor::all, j)) * v;
  }

  return ret;
}
/// @}

/// Outer product of two vectors
template <typename derived_1_t, typename derived_2_t>
DETRAY_HOST_DEVICE constexpr auto outer_product(
    const Fastor::AbstractTensor<derived_1_t, 1> &a,
    const Fastor::AbstractTensor<derived_2_t, 1> &b) {
  return Fastor::outer(a.self(), b.self());
}

/// @returns the determinant of @param m
template <std::size_t N, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr scalar_t determinant(
    const matrix_type<scalar_t, N, N> &m) {
  return Fastor::determinant(m);
}

/// @returns the inverse of @param m
template <std::size_t ROWS, std::size_t COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS> inverse(
    const matrix_type<scalar_t, ROWS, COLS> &m) {
  return Fastor::inverse(m);
}

}  // namespace detray::algebra::fastor::math
