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
  Intersection3D::Status status1 =
      (qe.first * qe.first < s_onSurfaceTolerance * s_onSurfaceTolerance)
          ? Intersection3D::Status::onSurface
          : Intersection3D::Status::reachable;

  if (bcheck and not isOnSurface(gctx, solution1, direction, bcheck)) {
    status1 = Intersection3D::Status::missed;
  }

  // Check the validity of the second solution
  Vector3D solution2 = position + qe.first * direction;
  Intersection3D::Status status2 =
      (qe.second * qe.second < s_onSurfaceTolerance * s_onSurfaceTolerance)
          ? Intersection3D::Status::onSurface
          : Intersection3D::Status::reachable;
  if (bcheck and not isOnSurface(gctx, solution2, direction, bcheck)) {
    status2 = Intersection3D::Status::missed;
  }

  const auto& tf = transform(gctx);
  // Set the intersection
  Intersection3D first(tf * solution1, qe.first, status1);
  Intersection3D second(tf * solution2, qe.second, status2);
  SurfaceIntersection cIntersection(first, this);
  // Check one if its valid or neither is valid
  bool check1 = status1 != Intersection3D::Status::missed or
                (status1 == Intersection3D::Status::missed and
                 status2 == Intersection3D::Status::missed);
  // Check and (eventually) go with the first solution
  if ((check1 and qe.first * qe.first < qe.second * qe.second) or
      status2 == Intersection3D::Status::missed) {
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

inline AlignmentRowVector ConeSurface::alignmentToPathDerivative(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  // The global position
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // The vector between position and center
  const ActsRowVector<AlignmentScalar, 3> pcRowVec =
      (position - center(gctx)).transpose();
  // The rotation
  const auto& rotation = transform(gctx).rotation();
  // The local frame x/y/z axis
  const Vector3D localXAxis = rotation.col(0);
  const Vector3D localYAxis = rotation.col(1);
  const Vector3D localZAxis = rotation.col(2);
  // The local coordinates
  const Vector3D localPos = rotation.transpose() * position;
  const double dx = direction.dot(localXAxis);
  const double dy = direction.dot(localYAxis);
  const double dz = direction.dot(localZAxis);
  // The normalization factor
  const double tanAlpha2 = bounds().tanAlpha() * bounds().tanAlpha();
  const double norm = 1. / (1. - dz * dz * (1 + tanAlpha2));
  // The direction transpose
  const ActsRowVector<AlignmentScalar, 3> dirRowVec = direction.transpose();
  // The derivative of path w.r.t. the local axes
  // @note The following calculations assume that the intersection of the track
  // with the cone always satisfy: localPos.z()*tanAlpha =perp(localPos)
  const ActsRowVector<AlignmentScalar, 3> localXAxisToPath =
      -2.0 * norm * (dx * pcRowVec + localPos.x() * dirRowVec);
  const ActsRowVector<AlignmentScalar, 3> localYAxisToPath =
      -2.0 * norm * (dy * pcRowVec + localPos.y() * dirRowVec);
  const ActsRowVector<AlignmentScalar, 3> localZAxisToPath =
      2.0 * norm * tanAlpha2 * (dz * pcRowVec + localPos.z() * dirRowVec) -
      4.0 * norm * norm * (1 + tanAlpha2) *
          (dx * localPos.x() + dy * localPos.y() -
           dz * localPos.z() * tanAlpha2) *
          dz * dirRowVec;
  // Calculate the derivative of local frame axes w.r.t its rotation
  const auto [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);
  // Initialize the derivative of propagation path w.r.t. local frame
  // translation (origin) and rotation
  AlignmentRowVector alignToPath = AlignmentRowVector::Zero();
  alignToPath.segment<3>(eAlignmentCenter0) =
      2.0 * norm * (dx * localXAxis.transpose() + dy * localYAxis.transpose());
  alignToPath.segment<3>(eAlignmentRotation0) =
      localXAxisToPath * rotToLocalXAxis + localYAxisToPath * rotToLocalYAxis +
      localZAxisToPath * rotToLocalZAxis;

  return alignToPath;
}

inline LocalCartesianToBoundLocalMatrix
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
