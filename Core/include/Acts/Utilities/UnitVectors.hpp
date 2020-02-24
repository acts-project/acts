// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <utility>

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// Construct the first curvilinear unit vector `U` for the given direction.
///
/// @param direction is the input direction vector
/// @returns a normalized vector in the x-y plane orthogonal to the direction.
///
/// The special case of the direction vector pointing along the z-axis is
/// handled by forcing the unit vector to along the x-axis.
template <typename InputVector>
inline auto makeCurvilinearUnitU(
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
inline auto makeCurvilinearUnitVectors(
    const Eigen::MatrixBase<InputVector>& direction) {
  EIGEN_STATIC_ASSERT_FIXED_SIZE(InputVector);
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(InputVector);
  static_assert(3 <= InputVector::RowsAtCompileTime,
                "Direction vector must be at least three-dimensional.");

  using OutputVector = typename InputVector::PlainObject;

  std::pair<OutputVector, OutputVector> unitVectors;
  unitVectors.first = makeCurvilinearUnitU(direction);
  unitVectors.second = direction.cross(unitVectors.first);
  unitVectors.second.normalize();
  return unitVectors;
}

}  // namespace Acts
