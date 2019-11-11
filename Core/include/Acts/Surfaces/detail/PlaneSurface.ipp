// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlaneSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

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
                                           const Vector3D& momentum) const {
  /// We can ignore the global position here
  return 1. /
         std::abs(Surface::normal(gctx, position).dot(momentum.normalized()));
}

inline Intersection PlaneSurface::intersectionEstimate(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  // Minimize the call to transform()
  const auto& tMatrix = transform(gctx).matrix();
  const Vector3D pnormal = tMatrix.block<3, 1>(0, 2).transpose();
  const Vector3D pcenter = tMatrix.block<3, 1>(0, 3).transpose();
  // return solution and path
  Vector3D solution(0., 0., 0.);
  double path = std::numeric_limits<double>::infinity();
  // Lemma : the solver -> should catch current values
  auto solve = [&solution, &path, &pnormal, &pcenter](
                   const Vector3D& lpos,
                   const Vector3D& ldir) -> Intersection::Status {
    double denom = ldir.dot(pnormal);
    if (denom != 0.0) {
      path = (pnormal.dot((pcenter - lpos))) / (denom);
      solution = (lpos + path * ldir);
      // Is valid hence either on surface or reachable
      return (path * path < s_onSurfaceTolerance * s_onSurfaceTolerance)
                 ? Intersection::Status::onSurface
                 : Intersection::Status::reachable;
    }
    return Intersection::Status::unreachable;
  };
  Intersection::Status status = solve(position, direction);
  // Evaluate (if necessary in terms of boundaries)
  // @todo: speed up isOnSurface - we know that it is on surface
  //  all we need is to check if it's inside bounds
  if (status != Intersection::Status::unreachable and bcheck and
      not isOnSurface(gctx, solution, direction, bcheck)) {
    status = Intersection::Status::missed;
  }
  return Intersection(solution, path, status);
}
