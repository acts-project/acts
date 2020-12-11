// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline Vector3 CylinderSurface::rotSymmetryAxis(
    const GeometryContext& gctx) const {
  // fast access via tranform matrix (and not rotation())
  return transform(gctx).matrix().block<3, 1>(0, 2);
}

inline detail::RealQuadraticEquation CylinderSurface::intersectionSolver(
    const Transform3& transform, const Vector3& position,
    const Vector3& direction) const {
  // Solve for radius R
  double R = bounds().get(CylinderBounds::eR);

  // Get the transformation matrtix
  const auto& tMatrix = transform.matrix();
  Vector3 caxis = tMatrix.block<3, 1>(0, 2).transpose();
  Vector3 ccenter = tMatrix.block<3, 1>(0, 3).transpose();

  // Check documentation for explanation
  Vector3 pc = position - ccenter;
  Vector3 pcXcd = pc.cross(caxis);
  Vector3 ldXcd = direction.cross(caxis);
  double a = ldXcd.dot(ldXcd);
  double b = 2. * (ldXcd.dot(pcXcd));
  double c = pcXcd.dot(pcXcd) - (R * R);
  // And solve the qaudratic equation
  return detail::RealQuadraticEquation(a, b, c);
}

inline SurfaceIntersection CylinderSurface::intersect(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const BoundaryCheck& bcheck) const {
  const auto& gctxTransform = transform(gctx);

  // Solve the quadratic equation
  auto qe = intersectionSolver(gctxTransform, position, direction);

  // If no valid solution return a non-valid surfaceIntersection
  if (qe.solutions == 0) {
    return SurfaceIntersection();
  }

  // Check the validity of the first solution
  Vector3 solution1 = position + qe.first * direction;
  Intersection3D::Status status1 =
      qe.first * qe.first < s_onSurfaceTolerance * s_onSurfaceTolerance
          ? Intersection3D::Status::onSurface
          : Intersection3D::Status::reachable;

  // Helper method for boundary check
  auto boundaryCheck =
      [&](const Vector3& solution,
          Intersection3D::Status status) -> Intersection3D::Status {
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
      const Vector3 vecLocal(solution - tMatrix.block<3, 1>(0, 3));
      double cZ = vecLocal.dot(tMatrix.block<3, 1>(0, 2));
      double tolerance = s_onSurfaceTolerance + bcheck.tolerance()[eBoundLoc1];
      double hZ = cBounds.get(CylinderBounds::eHalfLengthZ) + tolerance;
      return (cZ * cZ < hZ * hZ) ? status : Intersection3D::Status::missed;
    }
    return (isOnSurface(gctx, solution, direction, bcheck)
                ? status
                : Intersection3D::Status::missed);
  };
  // Check first solution for boundary compatiblity
  status1 = boundaryCheck(solution1, status1);
  // Set the intersection
  Intersection3D first(solution1, qe.first, status1);
  SurfaceIntersection cIntersection(first, this);
  if (qe.solutions == 1) {
    return cIntersection;
  }
  // Check the validity of the second solution
  Vector3 solution2 = position + qe.second * direction;
  Intersection3D::Status status2 =
      qe.second * qe.second < s_onSurfaceTolerance * s_onSurfaceTolerance
          ? Intersection3D::Status::onSurface
          : Intersection3D::Status::reachable;
  // Check first solution for boundary compatiblity
  status2 = boundaryCheck(solution2, status2);
  Intersection3D second(solution2, qe.second, status2);
  // Check one if its valid or neither is valid
  bool check1 = status1 != Intersection3D::Status::missed or
                (status1 == Intersection3D::Status::missed and
                 status2 == Intersection3D::Status::missed);
  // Check and (eventually) go with the first solution
  if ((check1 and qe.first * qe.first < qe.second * qe.second) or
      status2 == Intersection3D::Status::missed) {
    // And add the alternative
    cIntersection.alternative = second;
  } else {
    // And add the alternative
    cIntersection.alternative = first;
    cIntersection.intersection = second;
  }
  return cIntersection;
}

inline AlignmentToPathMatrix CylinderSurface::alignmentToPathDerivative(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  // The global position
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // The vector between position and center
  const auto pcRowVec = (position - center(gctx)).transpose().eval();
  // The rotation
  const auto& rotation = transform(gctx).rotation();
  // The local frame x/y/z axis
  const auto& localXAxis = rotation.col(0);
  const auto& localYAxis = rotation.col(1);
  const auto& localZAxis = rotation.col(2);
  // The local coordinates
  const auto localPos = (rotation.transpose() * position).eval();
  const auto dx = direction.dot(localXAxis);
  const auto dy = direction.dot(localYAxis);
  const auto dz = direction.dot(localZAxis);
  // The normalization factor
  const auto norm = 1 / (1 - dz * dz);
  // The direction transpose
  const auto& dirRowVec = direction.transpose();
  // The derivative of path w.r.t. the local axes
  // @note The following calculations assume that the intersection of the track
  // with the cylinder always satisfy: perp(localPos) = R
  const auto localXAxisToPath =
      (-2 * norm * (dx * pcRowVec + localPos.x() * dirRowVec)).eval();
  const auto localYAxisToPath =
      (-2 * norm * (dy * pcRowVec + localPos.y() * dirRowVec)).eval();
  const auto localZAxisToPath =
      (-4 * norm * norm * (dx * localPos.x() + dy * localPos.y()) * dz *
       dirRowVec)
          .eval();
  // Calculate the derivative of local frame axes w.r.t its rotation
  const auto [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);
  // Initialize the derivative of propagation path w.r.t. local frame
  // translation (origin) and rotation
  AlignmentToPathMatrix alignToPath = AlignmentToPathMatrix::Zero();
  alignToPath.segment<3>(eAlignmentCenter0) =
      2 * norm * (dx * localXAxis.transpose() + dy * localYAxis.transpose());
  alignToPath.segment<3>(eAlignmentRotation0) =
      localXAxisToPath * rotToLocalXAxis + localYAxisToPath * rotToLocalYAxis +
      localZAxisToPath * rotToLocalZAxis;

  return alignToPath;
}

inline ActsMatrix<2, 3> CylinderSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& gctx, const Vector3& position) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  // The local frame transform
  const auto& sTransform = transform(gctx);
  // calculate the transformation to local coorinates
  const Vector3 localPos = sTransform.inverse() * position;
  const double lr = perp(localPos);
  const double lphi = phi(localPos);
  const double lcphi = std::cos(lphi);
  const double lsphi = std::sin(lphi);
  // Solve for radius R
  double R = bounds().get(CylinderBounds::eR);
  ActsMatrix<2, 3> loc3DToLocBound = ActsMatrix<2, 3>::Zero();
  loc3DToLocBound << -R * lsphi / lr, R * lcphi / lr, 0, 0, 0, 1;

  return loc3DToLocBound;
}
