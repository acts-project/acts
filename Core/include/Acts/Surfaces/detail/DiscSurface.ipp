// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// DiscSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline const Vector2D DiscSurface::localPolarToCartesian(
    const Vector2D& lpolar) const {
  return Vector2D(lpolar[eLOC_R] * cos(lpolar[eLOC_PHI]),
                  lpolar[eLOC_R] * sin(lpolar[eLOC_PHI]));
}

inline const Vector2D DiscSurface::localCartesianToPolar(
    const Vector2D& lcart) const {
  return Vector2D(
      sqrt(lcart[eLOC_X] * lcart[eLOC_X] + lcart[eLOC_Y] * lcart[eLOC_Y]),
      atan2(lcart[eLOC_Y], lcart[eLOC_X]));
}

inline void DiscSurface::initJacobianToGlobal(const GeometryContext& gctx,
                                              BoundToFreeMatrix& jacobian,
                                              const Vector3D& position,
                                              const Vector3D& direction,
                                              const BoundVector& pars) const {
  // The trigonometry required to convert the direction to spherical
  // coordinates and then compute the sines and cosines again can be
  // surprisingly expensive from a performance point of view.
  //
  // Here, we can avoid it because the direction is by definition a unit
  // vector, with the following coordinate conversions...
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  const double z = direction(2);  // == cos(theta)

  // ...which we can invert to directly get the sines and cosines:
  const double cos_theta = z;
  const double sin_theta = sqrt(x * x + y * y);
  const double inv_sin_theta = 1. / sin_theta;
  const double cos_phi = x * inv_sin_theta;
  const double sin_phi = y * inv_sin_theta;
  // retrieve the reference frame
  const auto rframe = referenceFrame(gctx, position, direction);

  // special polar coordinates for the Disc
  double lrad = pars[eLOC_0];
  double lphi = pars[eLOC_1];
  double lcos_phi = cos(lphi);
  double lsin_phi = sin(lphi);
  // the local error components - rotated from reference frame
  jacobian.block<3, 1>(0, eLOC_0) =
      lcos_phi * rframe.block<3, 1>(0, 0) + lsin_phi * rframe.block<3, 1>(0, 1);
  jacobian.block<3, 1>(0, eLOC_1) =
      lrad * (lcos_phi * rframe.block<3, 1>(0, 1) -
              lsin_phi * rframe.block<3, 1>(0, 0));
  // the time component
  jacobian(3, eT) = 1;
  // the momentum components
  jacobian(4, ePHI) = (-sin_theta) * sin_phi;
  jacobian(4, eTHETA) = cos_theta * cos_phi;
  jacobian(5, ePHI) = sin_theta * cos_phi;
  jacobian(5, eTHETA) = cos_theta * sin_phi;
  jacobian(6, eTHETA) = (-sin_theta);
  jacobian(7, eQOP) = 1;
}

inline const RotationMatrix3D DiscSurface::initJacobianToLocal(
    const GeometryContext& gctx, FreeToBoundMatrix& jacobian,
    const Vector3D& position, const Vector3D& direction) const {
  using VectorHelpers::perp;
  using VectorHelpers::phi;
  // Optimized trigonometry on the propagation direction
  const double x = direction(0);  // == cos(phi) * sin(theta)
  const double y = direction(1);  // == sin(phi) * sin(theta)
  // component expressions - global
  const double inv_sin_theta_2 = 1. / (x * x + y * y);
  const double cos_phi_over_sin_theta = x * inv_sin_theta_2;
  const double sin_phi_over_sin_theta = y * inv_sin_theta_2;
  const double inv_sin_theta = sqrt(inv_sin_theta_2);
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
  jacobian.block<1, 3>(0, 0) = lcphi * lx + lsphi * ly;
  jacobian.block<1, 3>(1, 0) = (lcphi * ly - lsphi * lx) / lr;
  // Time element
  jacobian(eT, 3) = 1;
  // Directional and momentum elements for reference frame surface
  jacobian(ePHI, 4) = -sin_phi_over_sin_theta;
  jacobian(ePHI, 5) = cos_phi_over_sin_theta;
  jacobian(eTHETA, 6) = -inv_sin_theta;
  jacobian(eQOP, 7) = 1;
  // return the transposed reference frame
  return rframeT;
}

inline Intersection DiscSurface::intersectionEstimate(
    const GeometryContext& gctx, const Vector3D& position,
    const Vector3D& direction, const BoundaryCheck& bcheck) const {
  // minimize the call to transform()
  const auto& tMatrix = transform(gctx).matrix();
  const Vector3D pnormal = tMatrix.block<3, 1>(0, 2).transpose();
  const Vector3D pcenter = tMatrix.block<3, 1>(0, 3).transpose();
  // return solution and path
  Vector3D solution(0., 0., 0.);
  double path = std::numeric_limits<double>::infinity();
  // Lemma : the solver ------ encapsulated
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
  // --------
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

inline const Vector3D DiscSurface::normal(const GeometryContext& gctx,
                                          const Vector2D& /*unused*/) const {
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline const Vector3D DiscSurface::binningPosition(
    const GeometryContext& gctx, BinningValue /*unused*/) const {
  return center(gctx);
}

inline double DiscSurface::pathCorrection(const GeometryContext& gctx,
                                          const Vector3D& position,
                                          const Vector3D& momentum) const {
  /// we can ignore the global position here
  return 1. /
         std::abs(Surface::normal(gctx, position).dot(momentum.normalized()));
}
