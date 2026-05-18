// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/SpacePointFormation2/PixelSpacePointBuilder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

namespace Acts {

Vector2 PixelSpacePointBuilder::computeVarianceZR(
    const GeometryContext& gctx, const Surface& surface,
    const Vector3& spacePoint, const SquareMatrix2& localCov) {
  // the space point requires only the variance of the transverse and
  // longitudinal position. reduce computations by transforming the
  // covariance directly from local to z/r.
  //
  // compute Jacobian from global coordinates to z/r
  //
  //       dz/dz = 1
  //           r = sqrt(x² + y²)
  //   dr/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
  //             = 2 * {x,y} / r
  //
  const double x = spacePoint.x();
  const double y = spacePoint.y();
  const double scale = 2 / fastHypot(x, y);
  Matrix<2, 3> jacXyzToZr = Matrix<2, 3>::Zero();
  jacXyzToZr(0, 2) = 1;
  jacXyzToZr(1, 0) = scale * x;
  jacXyzToZr(1, 1) = scale * y;

  // using invalid direction vector, as it is usually not needed by the surface
  SquareMatrix3 rotLocalToGlobal =
      surface.referenceFrame(gctx, spacePoint, Vector3::Zero());

  // compute Jacobian from local coordinates to z/r
  SquareMatrix2 jac = jacXyzToZr * rotLocalToGlobal.topLeftCorner<3, 2>();

  // compute z/r variance
  Vector2 result = (jac * localCov * jac.transpose()).diagonal();
  return result;
}

}  // namespace Acts
