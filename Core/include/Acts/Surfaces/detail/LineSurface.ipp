// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// LineSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline void
LineSurface::localToGlobal(const Vector2D& lpos,
                           const Vector3D& mom,
                           Vector3D&       gpos) const
{
  // get the vector perpendicular to the momentum and the straw axis
  Vector3D radiusAxisGlobal(lineDirection().cross(mom));
  Vector3D locZinGlobal(0., 0., lpos[eLOC_Z]);
  // apply a transform if needed
  if (m_transform || (m_associatedDetElement != nullptr)) {
    locZinGlobal = transform() * locZinGlobal;
  }
  // transform zPosition into global coordinates and
  // add eLOC_R * radiusAxis
  gpos = Vector3D(locZinGlobal + lpos[eLOC_R] * radiusAxisGlobal.normalized());
}

inline bool
LineSurface::globalToLocal(const Vector3D& gpos,
                           const Vector3D& mom,
                           Vector2D&       lpos) const
{
  using VectorHelpers::perp;
  // apply the transform when needed
  Vector3D loc3Dframe = (m_transform || (m_associatedDetElement != nullptr))
      ? (transform().inverse()) * gpos
      : gpos;
  // construct localPosition with sign*perp(candidate) and z.()
  lpos = Vector2D(perp(loc3Dframe), loc3Dframe.z());
  Vector3D decVec(gpos - center());
  // assign the right sign
  double sign = ((lineDirection().cross(mom)).dot(decVec) < 0.) ? -1. : 1.;
  lpos[eLOC_R] *= sign;
  return true;
}

inline std::string
LineSurface::name() const
{
  return "Acts::LineSurface";
}

inline const RotationMatrix3D
LineSurface::referenceFrame(const Vector3D& /*pos*/, const Vector3D& mom) const
{
  RotationMatrix3D mFrame;
  // construct the measurement frame
  const Vector3D& measY = lineDirection();
  Vector3D        measX(measY.cross(mom).normalized());
  Vector3D        measDepth(measX.cross(measY));
  // assign the columnes
  mFrame.col(0) = measX;
  mFrame.col(1) = measY;
  mFrame.col(2) = measDepth;
  // return the rotation matrix
  return mFrame;
}

inline double
LineSurface::pathCorrection(const Vector3D& /*pos*/,
                            const Vector3D& /*mom*/) const
{
  return 1.;
}

inline const Vector3D
    LineSurface::binningPosition(BinningValue /*bValue*/) const
{
  return center();
}

inline const Vector3D
LineSurface::normal(const Vector2D& /*lpos*/) const
{
  // the normal is conceptionally closest to the line direction
  return lineDirection();
}

inline const SurfaceBounds&
LineSurface::bounds() const
{
  if (m_bounds) {
    return (*m_bounds.get());
  }
  return s_noBounds;
}

inline Intersection
LineSurface::intersectionEstimate(const Vector3D&      gpos,
                                  const Vector3D&      gdir,
                                  NavigationDirection  navDir,
                                  const BoundaryCheck& bcheck,
                                  CorrFnc              correct) const
{
  // following nominclature found in header file and doxygen documentation
  // line one is the straight track
  Vector3D ma = gpos;
  Vector3D ea = gdir;
  // line two is the line surface
  const auto& tMatrix = transform().matrix();
  Vector3D    mb      = tMatrix.block<3, 1>(0, 3).transpose();
  Vector3D    eb      = tMatrix.block<3, 1>(0, 2).transpose();
  // now go ahead and solve for the closest approach
  Vector3D mab(mb - ma);
  double   eaTeb = ea.dot(eb);
  double   denom = 1 - eaTeb * eaTeb;
  // validity parameter
  bool valid = false;
  if (denom * denom > s_onSurfaceTolerance * s_onSurfaceTolerance) {
    double u = (mab.dot(ea) - mab.dot(eb) * eaTeb) / denom;
    // evaluate in terms of direction
    valid = (navDir * u >= 0);
    // evaluate validaty in terms of bounds
    Vector3D result = (ma + u * ea);
    // update if you have a correction
    if (correct && correct(ma, ea, u)) {
      // update everything that is in relation to ea
      eaTeb = ea.dot(eb);
      denom = 1 - eaTeb * eaTeb;
      if (denom * denom > s_onSurfaceTolerance * s_onSurfaceTolerance) {
        u      = (mab.dot(ea) - mab.dot(eb) * eaTeb) / denom;
        result = (ma + u * ea);
        // if you have specified a navigation direction, valid mean path > 0.
        valid = (navDir * u >= 0);
      } else {
        valid = false;
      }
    }
    // it just needs to be a insideBounds() check
    // @todo there should be a faster check possible
    valid = bcheck ? (valid && isOnSurface(result, gdir, bcheck)) : valid;
    // return the result with validity
    return Intersection(result, u, valid);
  }
  // return a false intersection
  return Intersection(gpos, std::numeric_limits<double>::max(), false);
}

inline const Vector3D
LineSurface::lineDirection() const
{
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform().matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline void LineSurface::initJacobianToGlobal(ActsMatrixD<7, 5>& jacobian,
                                              const Vector3D&       gpos,
                                              const Vector3D&       dir,
                                              const ActsVectorD<5>& pars) const
{
  // The trigonometry required to convert the direction to spherical
  // coordinates and then compute the sines and cosines again can be
  // surprisingly expensive from a performance point of view.
  //
  // Here, we can avoid it because the direction is by definition a unit
  // vector, with the following coordinate conversions...
  const double x = dir(0);  // == cos(phi) * sin(theta)
  const double y = dir(1);  // == sin(phi) * sin(theta)
  const double z = dir(2);  // == cos(theta)

  // ...which we can invert to directly get the sines and cosines:
  const double cos_theta     = z;
  const double sin_theta     = sqrt(x * x + y * y);
  const double inv_sin_theta = 1. / sin_theta;
  const double cos_phi       = x * inv_sin_theta;
  const double sin_phi       = y * inv_sin_theta;
  // retrieve the reference frame
  const auto rframe = referenceFrame(gpos, dir);
  // the local error components - given by the reference frame
  jacobian.topLeftCorner<3, 2>() = rframe.topLeftCorner<3, 2>();
  // the momentum components
  jacobian(3, ePHI)   = (-sin_theta) * sin_phi;
  jacobian(3, eTHETA) = cos_theta * cos_phi;
  jacobian(4, ePHI)   = sin_theta * cos_phi;
  jacobian(4, eTHETA) = cos_theta * sin_phi;
  jacobian(5, eTHETA) = (-sin_theta);
  jacobian(6, eQOP)   = 1;

  // the projection of direction onto ref frame normal
  double ipdn = 1. / dir.dot(rframe.col(2));
  // build the cross product of d(D)/d(ePHI) components with y axis
  auto dDPhiY = rframe.block<3, 1>(0, 1).cross(jacobian.block<3, 1>(3, ePHI));
  // and the same for the d(D)/d(eTheta) components
  auto dDThetaY
      = rframe.block<3, 1>(0, 1).cross(jacobian.block<3, 1>(3, eTHETA));
  // and correct for the x axis components
  dDPhiY -= rframe.block<3, 1>(0, 0) * (rframe.block<3, 1>(0, 0).dot(dDPhiY));
  dDThetaY
      -= rframe.block<3, 1>(0, 0) * (rframe.block<3, 1>(0, 0).dot(dDThetaY));
  // set the jacobian components for global d/ phi/Theta
  jacobian.block<3, 1>(0, ePHI)   = dDPhiY * pars[eLOC_0] * ipdn;
  jacobian.block<3, 1>(0, eTHETA) = dDThetaY * pars[eLOC_0] * ipdn;
}

inline const ActsRowVectorD<5>
LineSurface::derivativeFactors(const Vector3D&         pos,
                               const Vector3D&         dir,
                               const RotationMatrix3D& rft,
                               const ActsMatrixD<7, 5>& jac) const
{
  // the vector between position and center
  ActsRowVectorD<3> pc = (pos - center()).transpose();
  // the longitudinal component vector (alogn local z)
  ActsRowVectorD<3> locz = rft.block<1, 3>(1, 0);
  // build the norm vector comonent by subtracting the longitudinal one
  double            long_c   = locz * dir;
  ActsRowVectorD<3> norm_vec = dir.transpose() - long_c * locz;
  // calculate the s factors for the dependency on X
  const ActsRowVectorD<5> s_vec = norm_vec * jac.topLeftCorner<3, 5>();
  // calculate the d factors for the dependency on Tx
  const ActsRowVectorD<5> d_vec = locz * jac.block<3, 5>(3, 0);
  // normalisation of normal & longitudinal components
  double norm = 1. / (1. - long_c * long_c);
  // create a matrix representation
  ActsMatrixD<3, 5> long_mat = ActsMatrixD<3, 5>::Zero();
  long_mat.colwise() += locz.transpose();
  // build the combined normal & longitudinal components
  return (
      norm
      * (s_vec - pc * (long_mat * d_vec.asDiagonal() - jac.block<3, 5>(3, 0))));
}
