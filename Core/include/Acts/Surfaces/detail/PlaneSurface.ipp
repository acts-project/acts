// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline Vector3 PlaneSurface::normal(const GeometryContext& gctx,
                                    const Vector2& /*lpos*/) const {
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline Vector3 PlaneSurface::binningPosition(const GeometryContext& gctx,
                                             BinningValue /*bValue*/) const {
  return center(gctx);
}

inline double PlaneSurface::pathCorrection(const GeometryContext& gctx,
                                           const Vector3& position,
                                           const Vector3& direction) const {
  // We can ignore the global position here
  return 1. / std::abs(Surface::normal(gctx, position).dot(direction));
}

inline SurfaceIntersection PlaneSurface::intersect(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryCheck& bcheck) const {
  // Get the contextual transform
  const auto& gctxTransform = transform(gctx);
  // Use the intersection helper for planar surfaces
  auto intersection =
      PlanarHelper::intersect(gctxTransform, position, direction);
  // Evaluate boundary check if requested (and reachable)
  if (intersection.status != Intersection3D::Status::unreachable and bcheck) {
    // Built-in local to global for speed reasons
    const auto& tMatrix = gctxTransform.matrix();
    // Create the reference vector in local
    const Vector3 vecLocal(intersection.position - tMatrix.block<3, 1>(0, 3));
    if (not insideBounds(tMatrix.block<3, 2>(0, 0).transpose() * vecLocal,
                         bcheck)) {
      intersection.status = Intersection3D::Status::missed;
    }
  }
  return {intersection, this};
}

inline ActsMatrix<2, 3> PlaneSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& /*unused*/, const Vector3& /*unused*/) const {
  const ActsMatrix<2, 3> loc3DToLocBound = ActsMatrix<2, 3>::Identity();
  return loc3DToLocBound;
}
