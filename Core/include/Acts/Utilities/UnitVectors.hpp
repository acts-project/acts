// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cmath>
#include <limits>
#include <utility>

namespace Acts {

/// Construct a normalized direction vector from phi angle and pseudorapidity.
///
/// @param phi is the direction angle in the x-y plane.
/// @param eta is the pseudorapidity towards the z-axis.
/// @return A normalized 3D direction vector constructed from phi and eta
///
/// @note The input arguments intentionally use the same template type so that
///       a compile error occurs if inconsistent input types are used. Avoids
///       unexpected implicit type conversions and forces the user to
///       explicitly cast mismatched input types.
template <typename T>
inline Eigen::Matrix<T, 3, 1> makeDirectionFromPhiEta(T phi, T eta) {
  const auto coshEtaInv = 1 / std::cosh(eta);
  return {
      std::cos(phi) * coshEtaInv,
      std::sin(phi) * coshEtaInv,
      std::tanh(eta),
  };
}

/// Construct a normalized direction vector from phi and theta angle.
///
/// @param phi is the direction angle in radian in the x-y plane.
/// @param theta is the polar angle in radian towards the z-axis.
/// @return A normalized 3D direction vector constructed from phi and theta
///
/// @note The input arguments intentionally use the same template type so that
///       a compile error occurs if inconsistent input types are used. Avoids
///       unexpected implicit type conversions and forces the user to
///       explicitly cast mismatched input types.
template <typename T>
inline Eigen::Matrix<T, 3, 1> makeDirectionFromPhiTheta(T phi, T theta) {
  const T sinTheta{std::sin(theta)};
  return {
      std::cos(phi) * sinTheta,
      std::sin(phi) * sinTheta,
      std::cos(theta),
  };
}
/// @brief Construct a normalized direction vector from the tangents of the
///        x-axis to the z-axis and of the y-axis to the z-axis
///
/// @param tanAlpha Tangent of the x-axis to the z-axis
/// @param tanBeta Tangent of the y-axis to the z-axis
/// @return A normalized 3D direction vector constructed from the axis tangents
template <typename T>
inline Eigen::Matrix<T, 3, 1> makeDirectionFromAxisTangents(T tanAlpha,
                                                            T tanBeta) {
  return Eigen::Matrix<T, 3, 1>{tanAlpha, tanBeta, 1}.normalized();
}

/// Construct a phi and theta angle from a direction vector.
///
/// @param unitDir 3D vector indicating a direction
/// @return A 2D vector containing phi and theta angles [phi, theta]
///
template <typename T>
inline Eigen::Matrix<T, 2, 1> makePhiThetaFromDirection(
    const Eigen::Matrix<T, 3, 1>& unitDir) {
  T phi = std::atan2(unitDir[1], unitDir[0]);
  T theta = std::acos(unitDir[2]);
  return {
      phi,
      theta,
  };
}

/// Construct the first curvilinear unit vector `U` for the given direction.
///
/// @param direction is the input direction vector
/// @returns a normalized vector in the x-y plane orthogonal to the direction.
///
/// The special case of the direction vector pointing along the z-axis is
/// handled by forcing the unit vector to along the x-axis.
template <typename InputVector>
inline auto createCurvilinearUnitU(
    const Eigen::MatrixBase<InputVector>& direction) {
  EIGEN_STATIC_ASSERT_FIXED_SIZE(InputVector);
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(InputVector);
  static_assert(3 <= InputVector::RowsAtCompileTime,
                "Direction vector must be at least three-dimensional.");

  using OutputVector = typename InputVector::PlainObject;
  using OutputScalar = typename InputVector::Scalar;

  OutputVector unitU = OutputVector::Zero();
  // explicit version of U = Z x T
  unitU[0] = -direction[1];
  unitU[1] = direction[0];
  const auto scale = unitU.template head<2>().norm();
  // if the absolute scale is tiny, the initial direction vector is aligned with
  // the z-axis. the ZxT product is ill-defined since any vector in the x-y
  // plane would be orthogonal to the direction. fix the U unit vector along the
  // x-axis to avoid this numerical instability.
  if (scale < (16 * std::numeric_limits<OutputScalar>::epsilon())) {
    unitU[0] = 1;
    unitU[1] = 0;
  } else {
    unitU.template head<2>() /= scale;
  }
  return unitU;
}

/// Construct the curvilinear unit vectors `U` and `V` for the given direction.
///
/// @param direction is the input direction vector
/// @returns normalized unit vectors `U` and `V` orthogonal to the direction.
///
/// With `T` the normalized input direction, the three vectors `U`, `V`, and
/// `T` form an orthonormal basis set, i.e. they satisfy
///
///     U x V = T
///     V x T = U
///     T x U = V
///
/// with the additional condition that `U` is located in the global x-y plane.
template <typename InputVector>
inline auto createCurvilinearUnitVectors(
    const Eigen::MatrixBase<InputVector>& direction) {
  EIGEN_STATIC_ASSERT_FIXED_SIZE(InputVector);
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(InputVector);
  static_assert(3 <= InputVector::RowsAtCompileTime,
                "Direction vector must be at least three-dimensional.");

  using OutputVector = typename InputVector::PlainObject;

  std::pair<OutputVector, OutputVector> unitVectors;
  unitVectors.first = createCurvilinearUnitU(direction);
  unitVectors.second = direction.cross(unitVectors.first);
  unitVectors.second.normalize();
  return unitVectors;
}

}  // namespace Acts
