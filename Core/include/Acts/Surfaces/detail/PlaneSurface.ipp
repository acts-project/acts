// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline const Vector3D PlaneSurface::normal(const GeometryContext& gctx,
                                           const Vector2D& /*lpos*/) const {
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline const Vector3D PlaneSurface::binningPosition(
    const GeometryContext& gctx, BinningValue /*bValue*/) const {
  return center(gctx);
}

inline double PlaneSurface::pathCorrection(const GeometryContext& gctx,
                                           const Vector3D& position,
                                           const Vector3D& direction) const {
  // We can ignore the global position here
  return 1. / std::abs(Surface::normal(gctx, position).dot(direction));
}

inline Intersection PlaneSurface::intersectionEstimate(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  // Get the contextual transform
  const auto& gctxTransform = transform(gctx);
  // Use the intersection helper for planar surfaces
  auto intersection =
      PlanarHelper::intersectionEstimate(gctxTransform, position, direction);
  // Evaluate boundary check if requested (and reachable)
  if (intersection.status != Intersection::Status::unreachable and bcheck) {
    // Built-in local to global for speed reasons
    const auto& tMatrix = gctxTransform.matrix();
    // Create the reference vector in local
    const Vector3D vecLocal(intersection.position - tMatrix.block<3, 1>(0, 3));
    if (not insideBounds(tMatrix.block<3, 2>(0, 0).transpose() * vecLocal,
                         bcheck)) {
      intersection.status = Intersection::Status::missed;
    }
  }
  return intersection;
}

inline const LocalCartesianToBoundLocalMatrix
PlaneSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& /*unused*/, const Vector3D& /*unused*/) const {
  const LocalCartesianToBoundLocalMatrix loc3DToLocBound =
      LocalCartesianToBoundLocalMatrix::Identity();
  return loc3DToLocBound;
}
