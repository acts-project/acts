// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

namespace Acts {

/// @brief Calculates the Jacobian for spherical to free
///        direction vector transformation
///
/// @note We use the direction vector as an input because
///       the trigonometric simplify that way
///
/// @param direction The normalised direction vector
///
/// @return The Jacobian d(dir_x, dir_y, dir_z) / d(phi, theta)
///
inline Matrix<3, 2> sphericalToFreeDirectionJacobian(const Vector3& direction) {
  auto [cosPhi, sinPhi, cosTheta, sinTheta] =
      VectorHelpers::evaluateTrigonomics(direction);

  // clang-format off
  Matrix<3, 2> jacobian;
  jacobian <<
    -direction.y(),  cosTheta * cosPhi,
     direction.x(),  cosTheta * sinPhi,
     0,             -sinTheta;
  // clang-format on

  return jacobian;
}

/// @brief Calculates the Jacobian for free to spherical
///        direction vector transformation
///
/// @note We use the direction vector as an input because
///       the trigonometric simplify that way
///
/// @param direction The normalised direction vector
///
/// @return The Jacobian d(phi, theta) / d(dir_x, dir_y, dir_z)
///
inline Matrix<2, 3> freeToSphericalDirectionJacobian(const Vector3& direction) {
  auto [cosPhi, sinPhi, cosTheta, sinTheta] =
      VectorHelpers::evaluateTrigonomics(direction);
  double invSinTheta = 1. / sinTheta;

  // clang-format off
  Matrix<2, 3> jacobian;
  jacobian <<
    -sinPhi * invSinTheta, cosPhi * invSinTheta, 0,
     cosPhi * cosTheta,    sinPhi * cosTheta,    -sinTheta;
  // clang-format on

  return jacobian;
}

/// @brief Calculates the Jacobian for spherical to free
///        momentum transformation
///
/// @note We use the direction vector as an input because
///       the trigonometric simplify that way
///
/// @param direction The normalised direction vector
/// @param qOverP The charge over momentum used in the free parametrization
/// @param momentum The absolute momentum matching @p qOverP for the
///        particle's charge hypothesis
///
/// @return The Jacobian d(p_x, p_y, p_z) / d(phi, theta, qOverP)
///
inline Matrix<3, 3> sphericalToFreeMomentumJacobian(const Vector3& direction,
                                                     double qOverP,
                                                     double momentum) {
  Matrix<3, 3> jacobian;
  jacobian.leftCols<2>() =
      momentum * sphericalToFreeDirectionJacobian(direction);
  jacobian.col(2) = -(momentum / qOverP) * direction;

  return jacobian;
}

/// @brief Calculates the Jacobian for free to spherical
///        momentum transformation
///
/// @param momentum The global, non-normalized momentum three-vector
/// @param charge The signed charge of the particle the momentum belongs to
///
/// @return The Jacobian d(phi, theta, qOverP) / d(p_x, p_y, p_z)
///
inline Matrix<3, 3> freeToSphericalMomentumJacobian(const Vector3& momentum,
                                                     double charge) {
  double p = momentum.norm();
  Vector3 direction = momentum / p;
  double qOverP = charge / p;

  Matrix<3, 3> jacobian;
  jacobian.topRows<2>() = freeToSphericalDirectionJacobian(direction) / p;
  jacobian.row(2) = -(qOverP / p) * direction.transpose();

  return jacobian;
}

}  // namespace Acts
