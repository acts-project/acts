// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/SpacePointFormation/PixelSpacePointBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

namespace Acts {

SquareMatrix2 PixelSpacePointBuilder::computeCovarianceZR(
    const RotationMatrix3& rotLocalToGlobal, const Vector3& spacePoint,
    const SquareMatrix2& localCov) {
  // the space point requires only the variance of the transverse and
  // longitudinal position. reduce computations by transforming the
  // covariance directly from local to z/r.
  //
  // compute Jacobian from global coordinates to z/r
  //
  //       dz/dz = 1
  //           r = sqrt(x² + y²)
  //   dr/d{x,y} = {x,y} / r
  //
  const double x = spacePoint.x();
  const double y = spacePoint.y();
  const double scale = 1 / fastHypot(x, y);
  Matrix<2, 3> jacXyzToZr = Matrix<2, 3>::Zero();
  jacXyzToZr(0, 2) = 1;
  jacXyzToZr(1, 0) = scale * x;
  jacXyzToZr(1, 1) = scale * y;

  // compute Jacobian from local coordinates to z/r
  const SquareMatrix2 jac = jacXyzToZr * rotLocalToGlobal.topLeftCorner<3, 2>();

  return jac * localCov * jac.transpose();
}

Vector2 PixelSpacePointBuilder::computeVarianceZR(
    const GeometryContext& gctx, const Surface& surface,
    const Vector3& spacePoint, const SquareMatrix2& localCov) {
  // using invalid direction vector, as it is usually not needed by the surface
  const RotationMatrix3 rotLocalToGlobal =
      surface.referenceFrame(gctx, spacePoint, Vector3::Zero());

  return computeCovarianceZR(rotLocalToGlobal, spacePoint, localCov).diagonal();
}

}  // namespace Acts
