// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/SurfaceBounds.hpp"

namespace Acts {

bool SurfaceBounds::inside(const Vector2& lposition,
                           const BoundaryTolerance& boundaryTolerance) const {
  using enum BoundaryTolerance::ToleranceMode;

  if (boundaryTolerance.isInfinite()) {
    return true;
  }

  BoundaryTolerance::ToleranceMode toleranceMode =
      boundaryTolerance.toleranceMode();
  bool strictlyInside = inside(lposition);

  if (toleranceMode == None) {
    return strictlyInside;
  }

  if (toleranceMode == Extend && strictlyInside) {
    return true;
  }

  SquareMatrix2 boundToCartesian = boundToCartesianJacobian(lposition);
  SquareMatrix2 metric = SquareMatrix2::Identity();
  if (boundaryTolerance.hasAbsoluteEuclidean()) {
    metric = boundToCartesian.transpose() * boundToCartesian;
  } else if (boundaryTolerance.hasChi2Bound()) {
    metric = boundToCartesian.transpose() *
             boundaryTolerance.asChi2Bound().weightMatrix() * boundToCartesian;
  } else if (boundaryTolerance.hasChi2Cartesian()) {
    metric = boundToCartesian.transpose() *
             boundaryTolerance.asChi2Cartesian().weightMatrix() *
             boundToCartesian;
  } else {
    throw std::runtime_error(
        "SurfaceBounds::inside: Unsupported boundary tolerance type.");
  }

  Vector2 closest = closestPoint(lposition, metric);
  Vector2 distance = closest - lposition;

  if (toleranceMode == Shrink) {
    return boundaryTolerance.isTolerated(distance, boundToCartesian) &&
           strictlyInside;
  }
  return boundaryTolerance.isTolerated(distance, boundToCartesian);
}

}  // namespace Acts
