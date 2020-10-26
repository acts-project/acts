// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline Vector2D DiscSurface::localPolarToCartesian(
    const Vector2D& lpolar) const {
  return Vector2D(lpolar[eBoundLoc0] * cos(lpolar[eBoundLoc1]),
                  lpolar[eBoundLoc0] * sin(lpolar[eBoundLoc1]));
}

inline Vector2D DiscSurface::localCartesianToPolar(
    const Vector2D& lcart) const {
  return Vector2D(sqrt(lcart[eBoundLoc0] * lcart[eBoundLoc0] +
                       lcart[eBoundLoc1] * lcart[eBoundLoc1]),
                  atan2(lcart[eBoundLoc1], lcart[eBoundLoc0]));
}

inline BoundToFreeMatrix DiscSurface::jacobianLocalToGlobal(
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
  // special polar coordinates for the Disc
  double lrad = boundParams[eBoundLoc0];
  double lphi = boundParams[eBoundLoc1];
  double lcos_phi = cos(lphi);
  double lsin_phi = sin(lphi);
  // retrieve the reference frame
  const auto rframe = referenceFrame(gctx, position, direction);
  // Initialize the jacobian from local to global
  BoundToFreeMatrix jacToGlobal = BoundToFreeMatrix::Zero();
  // the local error components - rotated from reference frame
  jacToGlobal.block<3, 1>(eFreePos0, eBoundLoc0) =
      lcos_phi * rframe.block<3, 1>(0, 0) + lsin_phi * rframe.block<3, 1>(0, 1);
  jacToGlobal.block<3, 1>(eFreePos0, eBoundLoc1) =
      lrad * (lcos_phi * rframe.block<3, 1>(0, 1) -
              lsin_phi * rframe.block<3, 1>(0, 0));
  // the time component
  jacToGlobal(eFreeTime, eBoundTime) = 1;
  // the momentum components
  jacToGlobal(eFreeDir0, eBoundPhi) = (-sin_theta) * sin_phi;
  jacToGlobal(eFreeDir0, eBoundTheta) = cos_theta * cos_phi;
  jacToGlobal(eFreeDir1, eBoundPhi) = sin_theta * cos_phi;
  jacToGlobal(eFreeDir1, eBoundTheta) = cos_theta * sin_phi;
  jacToGlobal(eFreeDir2, eBoundTheta) = (-sin_theta);
  jacToGlobal(eFreeQOverP, eBoundQOverP) = 1;
  return jacToGlobal;
}

inline FreeToBoundMatrix DiscSurface::jacobianGlobalToLocal(
    const GeometryContext& gctx, const FreeVector& parameters) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  // The global position
  const auto position = parameters.segment<3>(eFreePos0);
  // The direction
  const auto direction = parameters.segment<3>(eFreeDir0);
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // The measurement frame of the surface
  RotationMatrix3D rframeT =
      referenceFrame(gctx, position, direction).transpose();
  // calculate the transformation to local coorinates
  const Vector3D pos_loc = transform(gctx).inverse() * position;
  const double lr = perp(pos_loc);
  const double lphi = phi(pos_loc);
  const double lcphi = cos(lphi);
  const double lsphi = sin(lphi);
  // rotate into the polar coorindates
  auto lx = rframeT.block<1, 3>(0, 0);
  auto ly = rframeT.block<1, 3>(1, 0);
  // Initalize the jacobian from global to local
  FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
  // Local position component
  jacToLocal.block<1, 3>(eBoundLoc0, eFreePos0) = lcphi * lx + lsphi * ly;
  jacToLocal.block<1, 3>(eBoundLoc1, eFreePos0) =
      (lcphi * ly - lsphi * lx) / lr;
  // Time element
  jacToLocal(eBoundTime, eFreeTime) = 1;
  // Directional and momentum elements for reference frame surface
  jacToLocal(eBoundPhi, eFreeDir0) = -sinPhi * invSinTheta;
  jacToLocal(eBoundPhi, eFreeDir1) = cosPhi * invSinTheta;
  jacToLocal(eBoundTheta, eFreeDir0) = cosPhi * cosTheta;
  jacToLocal(eBoundTheta, eFreeDir1) = sinPhi * cosTheta;
  jacToLocal(eBoundTheta, eFreeDir2) = -sinTheta;
  jacToLocal(eBoundQOverP, eFreeQOverP) = 1;
  return jacToLocal;
}

inline SurfaceIntersection DiscSurface::intersect(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  // Get the contextual transform
  auto gctxTransform = transform(gctx);
  // Use the intersection helper for planar surfaces
  auto intersection =
      PlanarHelper::intersect(gctxTransform, position, direction);
  // Evaluate boundary check if requested (and reachable)
  if (intersection.status != Intersection3D::Status::unreachable and bcheck and
      m_bounds != nullptr) {
    // Built-in local to global for speed reasons
    const auto& tMatrix = gctxTransform.matrix();
    const Vector3D vecLocal(intersection.position - tMatrix.block<3, 1>(0, 3));
    const Vector2D lcartesian =
        tMatrix.block<3, 2>(0, 0).transpose() * vecLocal;
    if (bcheck.type() == BoundaryCheck::Type::eAbsolute and
        m_bounds->coversFullAzimuth()) {
      double tolerance = s_onSurfaceTolerance + bcheck.tolerance()[eBoundLoc0];
      if (not m_bounds->insideRadialBounds(VectorHelpers::perp(lcartesian),
                                           tolerance)) {
        intersection.status = Intersection3D::Status::missed;
      }
    } else if (not insideBounds(localCartesianToPolar(lcartesian), bcheck)) {
      intersection.status = Intersection3D::Status::missed;
    }
  }
  return {intersection, this};
}

inline LocalCartesianToBoundLocalMatrix
DiscSurface::localCartesianToBoundLocalDerivative(
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
  LocalCartesianToBoundLocalMatrix loc3DToLocBound =
      LocalCartesianToBoundLocalMatrix::Zero();
  loc3DToLocBound << lcphi, lsphi, 0, -lsphi / lr, lcphi / lr, 0;

  return loc3DToLocBound;
}

inline Vector3D DiscSurface::normal(const GeometryContext& gctx,
                                    const Vector2D& /*unused*/) const {
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline Vector3D DiscSurface::binningPosition(const GeometryContext& gctx,
                                             BinningValue bValue) const {
  if (bValue == binR) {
    double r = m_bounds->binningValueR();
    double phi = m_bounds->binningValuePhi();
    return Vector3D(r * cos(phi), r * sin(phi), center(gctx).z());
  }
  return center(gctx);
}

inline double DiscSurface::binningPositionValue(const GeometryContext& gctx,
                                                BinningValue bValue) const {
  // only modify binR
  if (bValue == binR) {
    return VectorHelpers::perp(center(gctx)) + m_bounds->binningValueR();
  }
  return GeometryObject::binningPositionValue(gctx, bValue);
}

inline double DiscSurface::pathCorrection(const GeometryContext& gctx,
                                          const Vector3D& position,
                                          const Vector3D& direction) const {
  /// we can ignore the global position here
  return 1. / std::abs(Surface::normal(gctx, position).dot(direction));
}
