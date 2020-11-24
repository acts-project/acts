// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline Vector3D LineSurface::localToGlobal(const GeometryContext& gctx,
                                           const Vector2D& lposition,
                                           const Vector3D& momentum) const {
  const auto& sTransform = transform(gctx);
  const auto& tMatrix = sTransform.matrix();
  Vector3D lineDirection(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));

  // get the vector perpendicular to the momentum and the straw axis
  Vector3D radiusAxisGlobal(lineDirection.cross(momentum));
  Vector3D locZinGlobal = sTransform * Vector3D(0., 0., lposition[eBoundLoc1]);
  // add eBoundLoc0 * radiusAxis
  return Vector3D(locZinGlobal +
                  lposition[eBoundLoc0] * radiusAxisGlobal.normalized());
}

inline Result<Vector2D> LineSurface::globalToLocal(const GeometryContext& gctx,
                                                   const Vector3D& position,
                                                   const Vector3D& momentum,
                                                   double /*tolerance*/) const {
  using VectorHelpers::perp;
  const auto& sTransform = transform(gctx);
  const auto& tMatrix = sTransform.matrix();
  Vector3D lineDirection(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
  // Bring the global position into the local frame
  Vector3D loc3Dframe = sTransform.inverse() * position;
  // construct localPosition with sign*perp(candidate) and z.()
  Vector2D lposition(perp(loc3Dframe), loc3Dframe.z());
  Vector3D sCenter(tMatrix(0, 3), tMatrix(1, 3), tMatrix(2, 3));
  Vector3D decVec(position - sCenter);
  // assign the right sign
  double sign = ((lineDirection.cross(momentum)).dot(decVec) < 0.) ? -1. : 1.;
  lposition[eBoundLoc0] *= sign;
  return Result<Vector2D>::success(lposition);
}

inline std::string LineSurface::name() const {
  return "Acts::LineSurface";
}

inline RotationMatrix3D LineSurface::referenceFrame(
    const GeometryContext& gctx, const Vector3D& /*unused*/,
    const Vector3D& momentum) const {
  RotationMatrix3D mFrame;
  const auto& tMatrix = transform(gctx).matrix();
  Vector3D measY(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
  Vector3D measX(measY.cross(momentum).normalized());
  Vector3D measDepth(measX.cross(measY));
  // assign the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  return mFrame;
}

inline double LineSurface::pathCorrection(const GeometryContext& /*unused*/,
                                          const Vector3D& /*pos*/,
                                          const Vector3D& /*mom*/) const {
  return 1.;
}

inline Vector3D LineSurface::binningPosition(const GeometryContext& gctx,
                                             BinningValue /*bValue*/) const {
  return center(gctx);
}

inline Vector3D LineSurface::normal(const GeometryContext& gctx,
                                    const Vector2D& /*lpos*/) const {
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline const SurfaceBounds& LineSurface::bounds() const {
  if (m_bounds) {
    return (*m_bounds.get());
  }
  return s_noBounds;
}

inline SurfaceIntersection LineSurface::intersect(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  // following nominclature found in header file and doxygen documentation
  // line one is the straight track
  const Vector3D& ma = position;
  const Vector3D& ea = direction;
  // line two is the line surface
  const auto& tMatrix = transform(gctx).matrix();
  Vector3D mb = tMatrix.block<3, 1>(0, 3).transpose();
  Vector3D eb = tMatrix.block<3, 1>(0, 2).transpose();
  // now go ahead and solve for the closest approach
  Vector3D mab(mb - ma);
  double eaTeb = ea.dot(eb);
  double denom = 1 - eaTeb * eaTeb;
  // validity parameter
  Intersection3D::Status status = Intersection3D::Status::unreachable;
  if (denom * denom > s_onSurfaceTolerance * s_onSurfaceTolerance) {
    double u = (mab.dot(ea) - mab.dot(eb) * eaTeb) / denom;
    // Check if we are on the surface already
    status = (u * u < s_onSurfaceTolerance * s_onSurfaceTolerance)
                 ? Intersection3D::Status::onSurface
                 : Intersection3D::Status::reachable;
    Vector3D result = (ma + u * ea);
    // Evaluate the boundary check if requested
    // m_bounds == nullptr prevents unecessary calulations for PerigeeSurface
    if (bcheck and m_bounds) {
      // At closest approach: check inside R or and inside Z
      const Vector3D vecLocal(result - mb);
      double cZ = vecLocal.dot(eb);
      double hZ =
          m_bounds->get(LineBounds::eHalfLengthZ) + s_onSurfaceTolerance;
      if ((cZ * cZ > hZ * hZ) or
          ((vecLocal - cZ * eb).norm() >
           m_bounds->get(LineBounds::eR) + s_onSurfaceTolerance)) {
        status = Intersection3D::Status::missed;
      }
    }
    return {Intersection3D(result, u, status), this};
  }
  // return a false intersection
  return {Intersection3D(position, std::numeric_limits<double>::max(), status),
          this};
}

inline BoundToFreeMatrix LineSurface::jacobianLocalToGlobal(
    const GeometryContext& gctx, const BoundVector& boundParams) const {
  // Transform from bound to free parameters
  FreeVector freeParams =
      detail::transformBoundToFreeParameters(*this, gctx, boundParams);
  // The global position
  const Vector3D position = freeParams.segment<3>(eFreePos0);
  // The direction
  const Vector3D direction = freeParams.segment<3>(eFreeDir0);
  // Get the sines and cosines directly
  const double cos_theta = std::cos(boundParams[eBoundTheta]);
  const double sin_theta = std::sin(boundParams[eBoundTheta]);
  const double cos_phi = std::cos(boundParams[eBoundPhi]);
  const double sin_phi = std::sin(boundParams[eBoundPhi]);
  // retrieve the reference frame
  const auto rframe = referenceFrame(gctx, position, direction);
  // Initialize the jacobian from local to global
  BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();
  // the local error components - given by the reference frame
  jacToGlobal.topLeftCorner<3, 2>() = rframe.topLeftCorner<3, 2>();
  // the time component
  jacToGlobal(eFreeTime, eBoundTime) = 1;
  // the momentum components
  jacToGlobal(eFreeDir0, eBoundPhi) = (-sin_theta) * sin_phi;
  jacToGlobal(eFreeDir0, eBoundTheta) = cos_theta * cos_phi;
  jacToGlobal(eFreeDir1, eBoundPhi) = sin_theta * cos_phi;
  jacToGlobal(eFreeDir1, eBoundTheta) = cos_theta * sin_phi;
  jacToGlobal(eFreeDir2, eBoundTheta) = (-sin_theta);
  jacToGlobal(eFreeQOverP, eBoundQOverP) = 1;

  // the projection of direction onto ref frame normal
  double ipdn = 1. / direction.dot(rframe.col(2));
  // build the cross product of d(D)/d(eBoundPhi) components with y axis
  auto dDPhiY = rframe.block<3, 1>(0, 1).cross(
      jacToGlobal.block<3, 1>(eFreeDir0, eBoundPhi));
  // and the same for the d(D)/d(eTheta) components
  auto dDThetaY = rframe.block<3, 1>(0, 1).cross(
      jacToGlobal.block<3, 1>(eFreeDir0, eBoundTheta));
  // and correct for the x axis components
  dDPhiY -= rframe.block<3, 1>(0, 0) * (rframe.block<3, 1>(0, 0).dot(dDPhiY));
  dDThetaY -=
      rframe.block<3, 1>(0, 0) * (rframe.block<3, 1>(0, 0).dot(dDThetaY));
  // set the jacobian components for global d/ phi/Theta
  jacToGlobal.block<3, 1>(eFreePos0, eBoundPhi) =
      dDPhiY * boundParams[eBoundLoc0] * ipdn;
  jacToGlobal.block<3, 1>(eFreePos0, eBoundTheta) =
      dDThetaY * boundParams[eBoundLoc0] * ipdn;
  return jacToGlobal;
}

inline FreeRowVector LineSurface::freeToPathDerivative(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  // The global posiiton
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // The vector between position and center
  const ActsRowVector<AlignmentScalar, 3> pcRowVec =
      (position - center(gctx)).transpose();
  // The rotation
  const auto& rotation = transform(gctx).rotation();
  // The local frame z axis
  const Vector3D localZAxis = rotation.col(2);
  // The local z coordinate
  const double pz = pcRowVec * localZAxis;
  // Cosine of angle between momentum direction and local frame z axis
  const double dz = localZAxis.dot(direction);
  const double norm = 1. / (1. - dz * dz);
  // Initialize the derivative of propagation path w.r.t. free parameter
  FreeRowVector freeToPath = FreeRowVector::Zero();
  // The derivative of path w.r.t. position
  freeToPath.segment<3>(eFreePos0) =
      norm * (dz * localZAxis.transpose() - direction.transpose());
  // The derivative of path w.r.t. direction
  freeToPath.segment<3>(eFreeDir0) =
      norm * (pz * localZAxis.transpose() - pcRowVec);

  return freeToPath;
}

inline AlignmentRowVector LineSurface::alignmentToPathDerivative(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  // The global posiiton
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // The vector between position and center
  const ActsRowVector<AlignmentScalar, 3> pcRowVec =
      (position - center(gctx)).transpose();
  // The rotation
  const auto& rotation = transform(gctx).rotation();
  // The local frame z axis
  const Vector3D localZAxis = rotation.col(2);
  // The local z coordinate
  const double pz = pcRowVec * localZAxis;
  // Cosine of angle between momentum direction and local frame z axis
  const double dz = localZAxis.dot(direction);
  const double norm = 1. / (1. - dz * dz);
  // Calculate the derivative of local frame axes w.r.t its rotation
  const auto [rotToLocalXAxis, rotToLocalYAxis, rotToLocalZAxis] =
      detail::rotationToLocalAxesDerivative(rotation);
  // Initialize the derivative of propagation path w.r.t. local frame
  // translation (origin) and rotation
  AlignmentRowVector alignToPath = AlignmentRowVector::Zero();
  alignToPath.segment<3>(eAlignmentCenter0) =
      norm * (direction.transpose() - dz * localZAxis.transpose());
  alignToPath.segment<3>(eAlignmentRotation0) =
      norm * (dz * pcRowVec + pz * direction.transpose()) * rotToLocalZAxis;

  return alignToPath;
}

inline LocalCartesianToBoundLocalMatrix
LineSurface::localCartesianToBoundLocalDerivative(
    const GeometryContext& gctx, const Vector3D& position) const {
  using VectorHelpers::phi;
  // The local frame transform
  const auto& sTransform = transform(gctx);
  // calculate the transformation to local coorinates
  const Vector3D localPos = sTransform.inverse() * position;
  const double lphi = phi(localPos);
  const double lcphi = std::cos(lphi);
  const double lsphi = std::sin(lphi);
  LocalCartesianToBoundLocalMatrix loc3DToLocBound =
      LocalCartesianToBoundLocalMatrix::Zero();
  loc3DToLocBound << lcphi, lsphi, 0, 0, 0, 1;

  return loc3DToLocBound;
}
