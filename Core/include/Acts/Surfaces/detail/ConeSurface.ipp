// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline detail::RealQuadraticEquation ConeSurface::intersectionSolver(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction) const {
  // Transform into the local frame
  Transform3D invTrans = transform(gctx).inverse();
  Vector3D point1 = invTrans * position;
  Vector3D dir1 = invTrans.linear() * direction;

  // See file header for the formula derivation
  double tan2Alpha = bounds().tanAlpha() * bounds().tanAlpha(),
         A = dir1.x() * dir1.x() + dir1.y() * dir1.y() -
             tan2Alpha * dir1.z() * dir1.z(),
         B = 2 * (dir1.x() * point1.x() + dir1.y() * point1.y() -
                  tan2Alpha * dir1.z() * point1.z()),
         C = point1.x() * point1.x() + point1.y() * point1.y() -
             tan2Alpha * point1.z() * point1.z();
  if (A == 0.) {
    A += 1e-16;  // avoid division by zero
  }

  return detail::RealQuadraticEquation(A, B, C);
}

inline Intersection ConeSurface::intersectionEstimate(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  // Solve the quadratic equation
  auto qe = intersectionSolver(gctx, position, direction);

  // If no valid solution return a non-valid intersection
  if (qe.solutions == 0) {
    return Intersection();
  }

  // Absolute smallest solution is taken by default
  double path =
      qe.first * qe.first < qe.second * qe.second ? qe.first : qe.second;
  Vector3D solution = position + path * direction;
  Intersection::Status status =
      (path * path < s_onSurfaceTolerance * s_onSurfaceTolerance)
          ? Intersection::Status::onSurface
          : Intersection::Status::reachable;

  // Boundary check necessary
  if (bcheck and not isOnSurface(gctx, solution, direction, bcheck)) {
    status = Intersection::Status::missed;
  }

  // Now return the solution
  return Intersection(transform(gctx) * solution, path, status);
}

inline SurfaceIntersection ConeSurface::intersect(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  // Solve the quadratic equation
  auto qe = intersectionSolver(gctx, position, direction);

  // If no valid solution return a non-valid surfaceIntersection
  if (qe.solutions == 0) {
    return SurfaceIntersection();
  }

  // Check the validity of the first solution
  Vector3D solution1 = position + qe.first * direction;
  Intersection::Status status1 =
      (qe.first * qe.first < s_onSurfaceTolerance * s_onSurfaceTolerance)
          ? Intersection::Status::onSurface
          : Intersection::Status::reachable;

  if (bcheck and not isOnSurface(gctx, solution1, direction, bcheck)) {
    status1 = Intersection::Status::missed;
  }

  // Check the validity of the second solution
  Vector3D solution2 = position + qe.first * direction;
  Intersection::Status status2 =
      (qe.second * qe.second < s_onSurfaceTolerance * s_onSurfaceTolerance)
          ? Intersection::Status::onSurface
          : Intersection::Status::reachable;
  if (bcheck and not isOnSurface(gctx, solution2, direction, bcheck)) {
    status2 = Intersection::Status::missed;
  }

  const auto& tf = transform(gctx);
  // Set the intersection
  Intersection first(tf * solution1, qe.first, status1);
  Intersection second(tf * solution2, qe.second, status2);
  SurfaceIntersection cIntersection(first, this);
  // Check one if its valid or neither is valid
  bool check1 = status1 != Intersection::Status::missed or
                (status1 == Intersection::Status::missed and
                 status2 == Intersection::Status::missed);
  // Check and (eventually) go with the first solution
  if ((check1 and qe.first * qe.first < qe.second * qe.second) or
      status2 == Intersection::Status::missed) {
    // And add the alternative
    if (qe.solutions > 1) {
      cIntersection.alternative = second;
    }
  } else {
    // And add the alternative
    if (qe.solutions > 1) {
      cIntersection.alternative = first;
    }
    cIntersection.intersection = second;
  }
  return cIntersection;
}

inline const LocalCartesianToBoundLocalMatrix
ConeSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& gctx, const Vector3D& position) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  // The local frame transform
  const auto& sTransform = transform(gctx);
  // calculate the transformation to local coorinates
  const Vector3D localPos = sTransform.inverse() * position;
  const double lr = perp(localPos);
  const double lphi = phi(localPos);
  const double lcphi = std::cos(lphi);
  const double lsphi = std::sin(lphi);
  // Solve for radius R
  const double R = localPos.z() * bounds().tanAlpha();
  LocalCartesianToBoundLocalMatrix loc3DToLocBound =
      LocalCartesianToBoundLocalMatrix::Zero();
  loc3DToLocBound << -R * lsphi / lr, R * lcphi / lr,
      lphi * bounds().tanAlpha(), 0, 0, 1;

  return loc3DToLocBound;
}
