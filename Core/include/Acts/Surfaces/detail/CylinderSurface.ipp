// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline const Vector3D CylinderSurface::rotSymmetryAxis(
    const GeometryContext& gctx) const {
  // fast access via tranform matrix (and not rotation())
  return transform(gctx).matrix().block<3, 1>(0, 2);
}

inline detail::RealQuadraticEquation CylinderSurface::intersectionSolver(
    const Transform3D& transform, const Vector3D& position,
    const Vector3D& direction) const {
  // Solve for radius R
  double R = bounds().get(CylinderBounds::eR);

  // Get the transformation matrtix
  const auto& tMatrix = transform.matrix();
  Vector3D caxis = tMatrix.block<3, 1>(0, 2).transpose();
  Vector3D ccenter = tMatrix.block<3, 1>(0, 3).transpose();

  // Check documentation for explanation
  Vector3D pc = position - ccenter;
  Vector3D pcXcd = pc.cross(caxis);
  Vector3D ldXcd = direction.cross(caxis);
  double a = ldXcd.dot(ldXcd);
  double b = 2. * (ldXcd.dot(pcXcd));
  double c = pcXcd.dot(pcXcd) - (R * R);
  // And solve the qaudratic equation
  return detail::RealQuadraticEquation(a, b, c);
}

inline Intersection CylinderSurface::intersectionEstimate(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  const auto& gctxTransform = transform(gctx);
  // Solve the quadratic equation
  auto qe = intersectionSolver(gctxTransform, position, direction);

  // If no valid solution return a non-valid intersection
  if (qe.solutions == 0) {
    return Intersection();
  }

  // Absolute smallest solution
  double path =
      qe.first * qe.first < qe.second * qe.second ? qe.first : qe.second;
  Vector3D solution = position + path * direction;
  Intersection::Status status =
      path * path < s_onSurfaceTolerance * s_onSurfaceTolerance
          ? Intersection::Status::onSurface
          : Intersection::Status::reachable;

  // Boundary check necessary
  if (bcheck and not isOnSurface(gctx, solution, direction, bcheck)) {
    status = Intersection::Status::missed;
  }

  // Now return the solution
  return Intersection(solution, path, status);
}

inline SurfaceIntersection CylinderSurface::intersect(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  const auto& gctxTransform = transform(gctx);

  // Solve the quadratic equation
  auto qe = intersectionSolver(gctxTransform, position, direction);

  // If no valid solution return a non-valid surfaceIntersection
  if (qe.solutions == 0) {
    return SurfaceIntersection();
  }

  // Check the validity of the first solution
  Vector3D solution1 = position + qe.first * direction;
  Intersection::Status status1 =
      qe.first * qe.first < s_onSurfaceTolerance * s_onSurfaceTolerance
          ? Intersection::Status::onSurface
          : Intersection::Status::reachable;

  // Helper method for boundary check
  auto boundaryCheck =
      [&](const Vector3D& solution,
          Intersection::Status status) -> Intersection::Status {
    // No check to be done, return current status
    if (!bcheck)
      return status;
    const auto& cBounds = bounds();
    if (cBounds.coversFullAzimuth() and
        bcheck.type() == BoundaryCheck::Type::eAbsolute) {
      // Project out the current Z value via local z axis
      // Built-in local to global for speed reasons
      const auto& tMatrix = gctxTransform.matrix();
      // Create the reference vector in local
      const Vector3D vecLocal(solution - tMatrix.block<3, 1>(0, 3));
      double cZ = vecLocal.dot(tMatrix.block<3, 1>(0, 2));
      double tolerance = s_onSurfaceTolerance + bcheck.tolerance()[eLOC_Z];
      double hZ = cBounds.get(CylinderBounds::eHalfLengthZ) + tolerance;
      return (cZ * cZ < hZ * hZ) ? status : Intersection::Status::missed;
    }
    return (isOnSurface(gctx, solution, direction, bcheck)
                ? status
                : Intersection::Status::missed);
  };
  // Check first solution for boundary compatiblity
  status1 = boundaryCheck(solution1, status1);
  // Set the intersection
  Intersection first(solution1, qe.first, status1);
  SurfaceIntersection cIntersection(first, this);
  if (qe.solutions == 1) {
    return cIntersection;
  }
  // Check the validity of the second solution
  Vector3D solution2 = position + qe.second * direction;
  Intersection::Status status2 =
      qe.second * qe.second < s_onSurfaceTolerance * s_onSurfaceTolerance
          ? Intersection::Status::onSurface
          : Intersection::Status::reachable;
  // Check first solution for boundary compatiblity
  status2 = boundaryCheck(solution2, status2);
  Intersection second(solution2, qe.second, status2);
  // Check one if its valid or neither is valid
  bool check1 = status1 != Intersection::Status::missed or
                (status1 == Intersection::Status::missed and
                 status2 == Intersection::Status::missed);
  // Check and (eventually) go with the first solution
  if ((check1 and qe.first * qe.first < qe.second * qe.second) or
      status2 == Intersection::Status::missed) {
    // And add the alternative
    cIntersection.alternative = second;
  } else {
    // And add the alternative
    cIntersection.alternative = first;
    cIntersection.intersection = second;
  }
  return cIntersection;
}

inline const LocalCartesianToBoundLocalMatrix
CylinderSurface::localCartesianToBoundLocalDerivative(
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
  double R = bounds().get(CylinderBounds::eR);
  LocalCartesianToBoundLocalMatrix loc3DToLocBound =
      LocalCartesianToBoundLocalMatrix::Zero();
  loc3DToLocBound << -R * lsphi / lr, R * lcphi / lr, 0, 0, 0, 1;

  return loc3DToLocBound;
}
