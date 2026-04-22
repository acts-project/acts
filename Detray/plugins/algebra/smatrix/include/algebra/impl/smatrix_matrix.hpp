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
#include "algebra/impl/smatrix_types.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// ROOT/Smatrix include(s).
#include <Math/SMatrix.h>

namespace detray::algebra::smatrix::math {

/// Create zero matrix
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr matrix_t zero() {
  return matrix_t();
}

/// Create identity matrix
template <concepts::matrix matrix_t>
DETRAY_HOST_DEVICE constexpr matrix_t identity() {
  return matrix_t(ROOT::Math::SMatrixIdentity());
}

/// Set input matrix as zero matrix
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr void set_zero(
    ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m) {
  m = zero<ROOT::Math::SMatrix<scalar_t, ROWS, COLS>>();
}

/// Set input matrix as identity matrix
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr void set_identity(
    ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m) {
  m = identity<ROOT::Math::SMatrix<scalar_t, ROWS, COLS>>();
}

/// Create transpose matrix
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS> transpose(
    const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m) {
  return ROOT::Math::Transpose(m);
}

/// Column-wise cross product
/// @{
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
  requires(ROWS == 3)
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, ROWS, COLS>
column_wise_cross(const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m,
                  const ROOT::Math::SVector<scalar_t, ROWS> &v) {
  matrix_type<scalar_t, ROWS, COLS> ret;
  for (unsigned int j = 0u; j < COLS; j++) {
    ROOT::Math::SVector<scalar_t, ROWS> col = ROOT::Math::Cross(
        m.template SubCol<ROOT::Math::SVector<scalar_t, ROWS>>(j, 0u), v);
    ret.Place_in_col(col, 0, j);
  }

  return ret;
}

template <typename OP, unsigned int ROWS, unsigned int COLS,
          concepts::scalar scalar_t,
          typename R = ROOT::Math::MatRepStd<scalar_t, ROWS, COLS>>
  requires(ROWS == 3)
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, ROWS, COLS>
column_wise_cross(const ROOT::Math::Expr<OP, scalar_t, ROWS, COLS, R> &m,
                  const ROOT::Math::SVector<scalar_t, ROWS> &v) {
  return column_wise_cross(static_cast<matrix_type<scalar_t, ROWS, COLS>>(m),
                           v);
}
/// @}

/// Column-wise product
/// @{
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, ROWS, COLS>
column_wise_multiply(const ROOT::Math::SMatrix<scalar_t, ROWS, COLS> &m,
                     const ROOT::Math::SVector<scalar_t, ROWS> &v) {
  matrix_type<scalar_t, ROWS, COLS> ret;
  for (unsigned int j = 0u; j < COLS; j++) {
    ROOT::Math::SVector<scalar_t, ROWS> col =
        m.template SubCol<ROOT::Math::SVector<scalar_t, ROWS>>(j, 0u) * v;
    ret.Place_in_col(col, 0, j);
  }

  return ret;
}

template <typename OP, unsigned int ROWS, unsigned int COLS,
          concepts::scalar scalar_t,
          typename R = ROOT::Math::MatRepStd<scalar_t, ROWS, COLS>>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, ROWS, COLS>
column_wise_multiply(const ROOT::Math::Expr<OP, scalar_t, ROWS, COLS, R> &m,
                     const ROOT::Math::SVector<scalar_t, ROWS> &v) {
  return column_wise_multiply(static_cast<matrix_type<scalar_t, ROWS, COLS>>(m),
                              v);
}
/// @}

/// Outer product of two vectors
template <unsigned int ROWS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, ROWS, ROWS> outer_product(
    const ROOT::Math::SVector<scalar_t, ROWS> &a,
    const ROOT::Math::SVector<scalar_t, ROWS> &b) {
  return ROOT::Math::TensorProd(a, b);
}

/// @returns the determinant of @param m
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr scalar_t determinant(
    const matrix_type<scalar_t, ROWS, COLS> &m) {
  scalar_t det;
  [[maybe_unused]] bool success = m.Det2(det);

  return det;
}

/// @returns the inverse of @param m
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, COLS, ROWS> inverse(
    const matrix_type<scalar_t, ROWS, COLS> &m) {
  int ifail = 0;
  return m.Inverse(ifail);
}

/// @returns the fatcor L in the decomposition of @param m
template <unsigned int ROWS, unsigned int COLS, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr matrix_type<scalar_t, ROWS, COLS>
cholesky_decomposition(const matrix_type<scalar_t, ROWS, COLS> &m) {
  matrix_type<scalar_t, ROWS, COLS> L;

  ROOT::Math::CholeskyDecomp<scalar_t, ROWS> decomp(m);
  [[maybe_unused]] const bool res = decomp.getL(L);
  assert(res);

  return L;
}

}  // namespace detray::algebra::smatrix::math
