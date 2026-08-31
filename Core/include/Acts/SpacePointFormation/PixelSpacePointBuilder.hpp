// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

namespace Acts {

class GeometryContext;
class Surface;

namespace PixelSpacePointBuilder {

/// @brief Transform a local covariance to z and r
///
/// Only for surfaces whose reference frame depends on the geometry context
/// alone, i.e. plane and disc, so that a caller can look it up once per module
/// with `Surface::referenceFrame(gctx)`.
///
/// @note Only the diagonal of the result feeds a space point. Prefer
/// @c computeVarianceZR, which skips the off-diagonal term and is inline.
///
/// @note A common cheaper shortcut is to split a single variance by the module
/// normal, `varZ = v * (1 - nz²)` and `varR = v * nz²`, dropping the second
/// local variance and the correlation.
///
/// @param rotLocalToGlobal The surface reference frame
/// @param spacePoint The global position
/// @param localCov The local covariance matrix
/// @return The (z, r) covariance
SquareMatrix2 computeCovarianceZR(const RotationMatrix3& rotLocalToGlobal,
                                  const Vector3& spacePoint,
                                  const SquareMatrix2& localCov);

/// @brief Transform a local covariance to the z and r variances
///
/// The diagonal of @c computeCovarianceZR, without forming the off-diagonal
/// term. Space point formation runs this once per measurement, so it is kept
/// inline and written out as two quadratic forms over the Jacobian rows.
///
/// See @c computeCovarianceZR for the requirement on the reference frame.
///
/// @param rotLocalToGlobal The surface reference frame
/// @param spacePoint The global position
/// @param localCov The local covariance matrix
/// @return The (z, r) variances
inline Vector2 computeVarianceZR(const RotationMatrix3& rotLocalToGlobal,
                                 const Vector3& spacePoint,
                                 const SquareMatrix2& localCov) {
  // the two local axes in global coordinates
  const double a0x = rotLocalToGlobal(0, 0);
  const double a0y = rotLocalToGlobal(1, 0);
  const double a0z = rotLocalToGlobal(2, 0);
  const double a1x = rotLocalToGlobal(0, 1);
  const double a1y = rotLocalToGlobal(1, 1);
  const double a1z = rotLocalToGlobal(2, 1);

  const double c00 = localCov(0, 0);
  const double c01 = localCov(0, 1) + localCov(1, 0);
  const double c11 = localCov(1, 1);

  // v^T localCov v for a two component Jacobian v
  const auto quadraticForm = [c00, c01, c11](double v0, double v1) {
    return c00 * v0 * v0 + c01 * v0 * v1 + c11 * v1 * v1;
  };

  const double x = spacePoint.x();
  const double y = spacePoint.y();

  // dz/dl is the z component of the local axis, dr/dl the radial one over r.
  // the radial Jacobian is left unnormalised, so its 1/r comes out of the
  // quadratic form as one division by r squared and no square root is needed
  const double varZ = quadraticForm(a0z, a1z);
  const double varR =
      quadraticForm(x * a0x + y * a0y, x * a1x + y * a1y) / (x * x + y * y);

  return {varZ, varR};
}

/// @brief Get z and r covariance from the local position and covariance
/// @param gctx The current geometry context object, e.g. alignment
/// @param surface The surface associated
/// @param spacePoint The global position
/// @param localCov The local covariance matrix
/// @return (z, r) components of the global covariance
/// @deprecated Use computeVarianceZR with a reference frame from
///             Surface::referenceFrame
[[deprecated(
    "Use computeVarianceZR with a reference frame from "
    "Surface::referenceFrame")]]
Vector2 computeVarianceZR(const GeometryContext& gctx, const Surface& surface,
                          const Vector3& spacePoint,
                          const SquareMatrix2& localCov);

}  // namespace PixelSpacePointBuilder

}  // namespace Acts
