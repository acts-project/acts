// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACE_SDISCSURFACE_H
#define ACTS_SURFACE_SDISCSURFACE_H 1

#include "ACTS/Surfaces/DiscBounds.hpp"
#include "ACTS/Surfaces/InfiniteBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Identifier.hpp"
#include "ACTS/Utilities/VariantDataFwd.hpp"

namespace Acts {

class DetectorElementBase;

/// @class DiscSurface
///
/// Class for a DiscSurface in the, it inherits from Surface.
///
/// The DiscSurface has a polar local coordinate system, with
/// (r,phi) describing the coordinates.
///
/// The surface transform positions the disc such, that the origin
/// is at r=0, independent of the provided DiscBounds. The z-axis
/// is the normal vector of the Disc, being perpendicular to the
/// radial direction.
///
/// The disc surface is the only surface type for which the
/// cavairance matrix is NOT given in the reference frame.
/// A conversion from polar to cartesian coordinates needs
/// to happen to transfer the local coordinates onto the
/// cartesian reference frame coordinates.
///
/// @image html DiscSurface.png
///
class DiscSurface : public Surface
{
public:
  /// Default Constructor is deleted
  DiscSurface() = delete;

  /// Constructor for Discs from Transform3D, \f$ r_{min}, r_{max} \f$
  ///
  /// @param htrans is transform that places the disc in the global 3D space
  /// (can be nullptr)
  /// @param rmin is the inner radius of the disc surface
  /// @param rmax is the outer radius of the disc surface
  /// @param hphisec is the opening angle of the disc surface and is optional
  ///        the default is a full disc
  DiscSurface(std::shared_ptr<const Transform3D> htrans,
              double                             rmin,
              double                             rmax,
              double                             hphisec = M_PI);

  /// Constructor for Discs from Transform3D, \f$ r_{min}, r_{max}, hx_{min},
  /// hx_{max} \f$
  /// This is n this case you have DiscTrapezoidalBounds
  ///
  /// @param htrans is transform that places the disc in the global 3D space
  /// (can be nullptr)
  /// @param minhalfx is the half length in x at minimal r
  /// @param maxhalfx is the half length in x at maximal r
  /// @param rmin is the inner radius of the disc surface
  /// @param rmax is the outer radius of the disc surface
  /// @param avephi is the position in phi (default is 0.)
  /// @param stereo is the optional stereo angle
  DiscSurface(std::shared_ptr<const Transform3D> htrans,
              double                             minhalfx,
              double                             maxhalfx,
              double                             rmin,
              double                             rmax,
              double                             avephi = 0.,
              double                             stereo = 0.);

  /// Constructor for Discs from Transform3D and shared DiscBounds
  ///
  /// @param htrans is the transform that positions the disc in the global 3D
  /// frame
  /// @param dbounds are the disc bounds describing the surface coverage
  DiscSurface(std::shared_ptr<const Transform3D> htrans,
              std::shared_ptr<const DiscBounds>  dbounds = nullptr);

  /// Constructor from detector element and identifier
  /// @note the surface only acts as a proxy of the detector element
  ///
  /// @param dbounds are the disc bounds associated to this surface, must not be
  /// nullptr
  /// @param detelement is the detector element that is represented by this
  /// surface
  /// @param identifier is the optional identifier in case one detector element
  /// owns more than 1 surface
  DiscSurface(std::shared_ptr<const DiscBounds> dbounds,
              const DetectorElementBase&        detelement,
              const Identifier&                 identifier = Identifier());

  /// Copy Constructor
  ///
  /// @param other is the source surface for the copy
  DiscSurface(const DiscSurface& other);

  /// Copy Constructor with shift
  ///
  /// @param other is the source sourface for the copy
  /// @param transf is the additional transform applied to the surface
  DiscSurface(const DiscSurface& other, const Transform3D& transf);

  DiscSurface(const variant_data& data);

  virtual ~DiscSurface();

  /// Assignement operator
  ///
  /// @param other is the source sourface for the assignment
  DiscSurface&
  operator=(const DiscSurface& other);

  /// Virtual constructor - shift can be given optionally
  ///
  /// @param shift the otional transform applied after cloning
  virtual DiscSurface*
  clone(const Transform3D* shift = nullptr) const final override;

  /// Return the surface type
  virtual SurfaceType
  type() const override;

  /// Normal vector
  ///
  /// @param lpos the local position where the normal is requested (ignored)
  /// @return is a normal vector
  const Vector3D
  normal(const Vector2D& lpos = s_origin2D) const final;

  /// The binning position is the position calcualted
  /// for a certain binning type
  ///
  /// @param bValue is the binning type to be used
  /// @return position that can beused for this binning
  virtual const Vector3D
  binningPosition(BinningValue bValue) const final;

  /// This method returns the bounds by reference
  const SurfaceBounds&
  bounds() const final override;

  /// This method returns true if the GlobalPosition is on the Surface for both,
  /// within or without check of whether the local position is inside boundaries
  /// or not
  ///
  /// @param gpos is the global position to be checked
  /// @param bcheck is the boundary check directive
  /// @return bollean that indicates if the position is on surface
  virtual bool
  isOnSurface(const Vector3D&      gpos,
              const BoundaryCheck& bcheck = true) const final override;

  /// Local to global transformation
  /// For planar surfaces the momentum is ignroed in the local to global
  /// transformation
  ///
  /// @param lpos local 2D posittion in specialized surface frame
  /// @param mom global 3D momentum representation (optionally ignored)
  /// @param gpos global 3D position to be filled (given by reference for method
  /// symmetry)
  ///
  /// @note the momentum is ignored for Disc surfaces in this calculateion
  virtual void
  localToGlobal(const Vector2D& lpos,
                const Vector3D& mom,
                Vector3D&       gpos) const final override;

  /// Global to local transformation
  /// @note the momentum is ignored for Disc surfaces in this calculateion
  ///
  /// @param gpos global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param mom global 3D momentum representation (optionally ignored)
  /// @param lpos local 2D position to be filled (given by reference for method
  /// symmetry)
  /// @return boolean indication if operation was successful (fail means global
  /// position was not on surface)
  virtual bool
  globalToLocal(const Vector3D& gpos,
                const Vector3D& mom,
                Vector2D&       lpos) const final override;

  /// Special method for DiscSurface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lpolar is a local position in polar coordinates
  /// @return values is local 2D position in cartesian coordinates  @todo check
  const Vector2D
  localPolarToCartesian(const Vector2D& lpolar) const;

  /// Special method for Disc surface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lcart is local 2D position in cartesian coordinates
  /// @return value is a local position in polar coordinates
  const Vector2D
  localCartesianToPolar(const Vector2D& lcart) const;

  /// Special method for DiscSurface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lpolar is a local position in polar coordinates
  /// @return values is local 2D position in cartesian coordinates
  const Vector2D
  localPolarToLocalCartesian(const Vector2D& lpolar) const;

  /// Special method for DiscSurface :  local<->global transformation when
  /// provided cartesian coordinates
  ///
  /// @param lcart is local 2D position in cartesian coordinates
  /// @return value is a global cartesian 3D position
  const Vector3D
  localCartesianToGlobal(const Vector2D& lcart) const;

  /// Special method for DiscSurface : global<->local from cartesian coordinates
  ///
  /// @param gpos is a global cartesian 3D position
  /// @param tol is the absoltue tolerance parameter
  /// @return value is a local polar
  const Vector2D
  globalToLocalCartesian(const Vector3D& gpos, double tol = 0.) const;

  /// Initialize the jacobian from local to global
  /// the surface knows best, hence the calculation is done here.
  /// The jacobian is assumed to be initialised, so only the
  /// relevant entries are filled
  ///
  /// @param jac is the jacobian to be initialized
  /// @param pos is the global position of the parameters
  /// @param dir is the direction at of the parameters
  /// @param pars is the paranmeters vector
  virtual void
      initJacobianToGlobal(ActsMatrixD<7, 5>& jac,
                           const Vector3D&       gpos,
                           const Vector3D&       dir,
                           const ActsVectorD<5>& pars) const final override;

  /// Initialize the jacobian from global to local
  /// the surface knows best, hence the calculation is done here.
  /// The jacobian is assumed to be initialised, so only the
  /// relevant entries are filled
  ///
  /// @param jac is the jacobian to be initialized
  /// @param pos is the global position of the parameters
  /// @param dir is the direction at of the parameters
  /// @param pars is the parameter vector
  ///
  /// @return the transposed reference frame (avoids recalculation)
  virtual const RotationMatrix3D
      initJacobianToLocal(ActsMatrixD<5, 7>& jac,
                          const Vector3D& gpos,
                          const Vector3D& dir) const final override;

  /// Path correction due to incident of the track
  ///
  /// @param gpos is the global position as a starting point
  /// @param mom is the global momentum at the starting point
  /// @return is the correction factor due to incident
  double
  pathCorrection(const Vector3D& gpos,
                 const Vector3D& mom) const final override;

  /// fast straight line intersection schema - standard: provides closest
  /// intersection and (signed) path length
  ///  forceDir is to provide the closest forward solution
  ///
  /// @param gpos is the global position as a starting point
  /// @param gdir is the global direction at the starting point
  ///        @note has to be normalized
  /// @param forceDir is a boolean forcing a solution along direction
  /// @param bcheck is the boundary check
  ///
  ///  <b>mathematical motivation:</b>
  ///
  ///  the equation of the plane is given by: <br>
  ///  @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
  ///  where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of
  /// the plane,
  ///  @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point on the plane and
  /// @f$ \vec x = (x,y,z) @f$ all possible points
  ///  on the plane.<br>
  ///  Given a line with:<br>
  ///  @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
  ///  the solution for @f$ u @f$ can be written:
  ///  @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
  ///  If the denominator is 0 then the line lies:
  ///  - either in the plane
  ///  - perpendicular to the normal of the plane
  ///
  /// @return is the intersection object
  virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      gdir,
                       bool                 forceDir = false,
                       const BoundaryCheck& bcheck
                       = false) const final override;

  /// Return properly formatted class name for screen output
  virtual std::string
  name() const override;

  virtual variant_data
  toVariantData() const override;

protected:
  std::shared_ptr<const DiscBounds> m_bounds;  ///< bounds (shared)
};

inline const Vector2D
DiscSurface::localPolarToCartesian(const Vector2D& lpolar) const
{
  return Vector2D(lpolar[Acts::eLOC_R] * cos(lpolar[Acts::eLOC_PHI]),
                  lpolar[Acts::eLOC_R] * sin(lpolar[Acts::eLOC_PHI]));
}

inline const Vector2D
DiscSurface::localCartesianToPolar(const Vector2D& lcart) const
{
  return Vector2D(sqrt(lcart[Acts::eLOC_X] * lcart[Acts::eLOC_X]
                       + lcart[Acts::eLOC_Y] * lcart[Acts::eLOC_Y]),
                  atan2(lcart[Acts::eLOC_Y], lcart[Acts::eLOC_X]));
}

inline void DiscSurface::initJacobianToGlobal(ActsMatrixD<7, 5>& jacobian,
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

  // special polar coordinates for the Disc
  double lrad     = pars[eLOC_0];
  double lphi     = pars[eLOC_1];
  double lcos_phi = cos(lphi);
  double lsin_phi = sin(lphi);
  // the local error components - rotated from reference frame
  jacobian.block<3, 1>(0, eLOC_0) = lcos_phi * rframe.block<3, 1>(0, 0)
      + lsin_phi * rframe.block<3, 1>(0, 1);
  jacobian.block<3, 1>(0, eLOC_1)
      = lrad * (lcos_phi * rframe.block<3, 1>(0, 1)
                - lsin_phi * rframe.block<3, 1>(0, 0));
  // the momentum components
  jacobian(3, ePHI)   = (-sin_theta) * sin_phi;
  jacobian(3, eTHETA) = cos_theta * cos_phi;
  jacobian(4, ePHI)   = sin_theta * cos_phi;
  jacobian(4, eTHETA) = cos_theta * sin_phi;
  jacobian(5, eTHETA) = (-sin_theta);
  jacobian(6, eQOP)   = 1;
}

inline const RotationMatrix3D
    DiscSurface::initJacobianToLocal(ActsMatrixD<5, 7>& jacobian,
                                     const Vector3D& gpos,
                                     const Vector3D& dir) const
{
  // Optimized trigonometry on the propagation direction
  const double x = dir(0);  // == cos(phi) * sin(theta)
  const double y = dir(1);  // == sin(phi) * sin(theta)
  // component expressions - global
  const double inv_sin_theta_2        = 1. / (x * x + y * y);
  const double cos_phi_over_sin_theta = x * inv_sin_theta_2;
  const double sin_phi_over_sin_theta = y * inv_sin_theta_2;
  const double inv_sin_theta          = sqrt(inv_sin_theta_2);
  // The measurement frame of the surface
  RotationMatrix3D rframeT = referenceFrame(gpos, dir).transpose();
  // calculate the transformation to local coorinates
  const Vector3D pos_loc = transform().inverse() * gpos;
  const double   lr      = pos_loc.perp();
  const double   lphi    = pos_loc.phi();
  const double   lcphi   = cos(lphi);
  const double   lsphi   = sin(lphi);
  // rotate into the polar coorindates
  auto lx = rframeT.block<1, 3>(0, 0);
  auto ly = rframeT.block<1, 3>(1, 0);
  jacobian.block<1, 3>(0, 0) = lcphi * lx + lsphi * ly;
  jacobian.block<1, 3>(1, 0) = (lcphi * ly - lsphi * lx) / lr;
  // Directional and momentum elements for reference frame surface
  jacobian(ePHI, 3)   = -sin_phi_over_sin_theta;
  jacobian(ePHI, 4)   = cos_phi_over_sin_theta;
  jacobian(eTHETA, 5) = -inv_sin_theta;
  jacobian(eQOP, 6)   = 1;
  // return the transposed reference frame
  return rframeT;
}

}  // end of namespace

#endif  // ACTS_SURFACES_DISCSURFACE_H
