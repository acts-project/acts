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
  // the Jacobian from global coordinates to z/r is
  //
  //       dz/dz = 1
  //           r = sqrt(x² + y²)
  //   dr/d{x,y} = {x,y} / r
  //
  // and is applied to the reference frame directly rather than multiplied out,
  // since three of its six entries are structurally zero.
  const double x = spacePoint.x();
  const double y = spacePoint.y();
  const double scale = 1 / fastHypot(x, y);

  // Jacobian from the two local coordinates to z/r: dz/dl is the z component
  // of the local axis, dr/dl its radial component
  SquareMatrix2 jac;
  jac(0, 0) = rotLocalToGlobal(2, 0);
  jac(0, 1) = rotLocalToGlobal(2, 1);
  jac(1, 0) = scale * (x * rotLocalToGlobal(0, 0) + y * rotLocalToGlobal(1, 0));
  jac(1, 1) = scale * (x * rotLocalToGlobal(0, 1) + y * rotLocalToGlobal(1, 1));

  return jac * localCov * jac.transpose();
}

Vector2 PixelSpacePointBuilder::computeVarianceZR(
    const GeometryContext& gctx, const Surface& surface,
    const Vector3& spacePoint, const SquareMatrix2& localCov) {
  // using invalid direction vector, as it is usually not needed by the surface
  const RotationMatrix3 rotLocalToGlobal =
      surface.referenceFrame(gctx, spacePoint, Vector3::Zero());

  return computeVarianceZR(rotLocalToGlobal, spacePoint, localCov);
}

}  // namespace Acts
