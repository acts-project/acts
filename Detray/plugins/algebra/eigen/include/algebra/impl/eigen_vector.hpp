// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
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

namespace detray::algebra::eigen::math {

/// This method retrieves phi from a vector @param v
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr auto phi(
    const Eigen::MatrixBase<derived_type> &v) {
  return algebra::math::atan2(v[1], v[0]);
}

/// This method retrieves the perpendicular magnitude of a vector @param v
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr auto perp(
    const Eigen::MatrixBase<derived_type> &v) {
  return algebra::math::sqrt(algebra::math::fma(v[0], v[0], v[1] * v[1]));
}

/// This method retrieves theta from a vector @param v
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr auto theta(
    const Eigen::MatrixBase<derived_type> &v) {
  return algebra::math::atan2(perp(v), v[2]);
}

/// This method retrieves the norm of a vector, no dimension restriction
///
/// @param v the input vector
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr auto norm(
    const Eigen::MatrixBase<derived_type> &v) {
  return v.norm();
}

/// This method retrieves the pseudo-rapidity from a vector or vector base with
/// rows >= 3
///
/// @param v the input vector
template <typename derived_type>
  requires(Eigen::MatrixBase<derived_type>::RowsAtCompileTime >= 3)
DETRAY_HOST_DEVICE constexpr auto eta(
    const Eigen::MatrixBase<derived_type> &v) noexcept {
  return algebra::math::atanh(v[2] / v.norm());
}

/// Get a normalized version of the input vector
///
/// @tparam derived_type is the matrix template
///
/// @param v the input vector
template <typename derived_type>
DETRAY_HOST_DEVICE constexpr auto normalize(
    const Eigen::MatrixBase<derived_type> &v) {
  return v.normalized();
}

/// Dot product between two input vectors
///
/// @tparam derived_type_lhs is the first matrix (expression) template
/// @tparam derived_type_rhs is the second matrix (expression) template
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return the scalar dot product value
template <typename derived_type_lhs, typename derived_type_rhs>
DETRAY_HOST_DEVICE constexpr auto dot(
    const Eigen::MatrixBase<derived_type_lhs> &a,
    const Eigen::MatrixBase<derived_type_rhs> &b) {
  return a.dot(b);
}

/// Cross product between two input vectors
///
/// @tparam derived_type_lhs is the first matrix (expression) template
/// @tparam derived_type_rhs is the second matrix (expression) template
///
/// @param a the first input vector
/// @param b the second input vector
///
/// @return a vector (expression) representing the cross product
template <typename derived_type_lhs, typename derived_type_rhs>
DETRAY_HOST_DEVICE constexpr auto cross(
    const Eigen::MatrixBase<derived_type_lhs> &a,
    const Eigen::MatrixBase<derived_type_rhs> &b) {
  return a.cross(b);
}

}  // namespace detray::algebra::eigen::math
