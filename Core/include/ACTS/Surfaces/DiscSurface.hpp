// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
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
/// is the normal vector of the Disc, being perpenticular to the
/// radial direction.
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
  DiscSurface(std::shared_ptr<Transform3D> htrans,
              double                       rmin,
              double                       rmax,
              double                       hphisec = M_PI);

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
  DiscSurface(std::shared_ptr<Transform3D> htrans,
              double                       minhalfx,
              double                       maxhalfx,
              double                       rmin,
              double                       rmax,
              double                       avephi = 0.,
              double                       stereo = 0.);

  /// Constructor for Discs from Transform3D and shared DiscBounds
  ///
  /// @param htrans is the transform that positions the disc in the global 3D
  /// frame
  /// @param dbounds are the disc bounds describing the surface coverage
  DiscSurface(std::shared_ptr<Transform3D>      htrans,
              std::shared_ptr<const DiscBounds> dbounds = nullptr);

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
  /// @param dsf is the source surface for the copy
  DiscSurface(const DiscSurface& dsf);

  /// Copy Constructor with shift
  ///
  /// @param dsf is the source sourface for the copy
  /// @param transf is the additional transform applied to the surface
  DiscSurface(const DiscSurface& dsf, const Transform3D& transf);

  /// Destructor
  virtual ~DiscSurface();

  /// Assignement operator
  ///
  /// @param dsf is the source sourface for the assignment
  DiscSurface&
  operator=(const DiscSurface& dsf);

  /// Virtual constructor - shift can be given optionally
  ///
  /// @param shift the otional transform applied after cloning
  virtual DiscSurface*
  clone(const Transform3D* shift = nullptr) const final override;

  /// Return the surface type
  virtual SurfaceType
  type() const override
  {
    return Surface::Disc;
  }

  /// Normal vector
  ///
  /// @param lpos the local position where the normal is requested (ignored)
  ///
  /// @return is a normal vector
  const Vector3D
  normal(const Vector2D& lpos = s_origin2D) const final;

  /// The binning position is the position calcualted
  /// for a certain binning type
  ///
  /// @param bValue is the binning type to be used
  ///
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
  ///
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
  ///
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
  ///
  /// @return values is local 2D position in carthesian coordinates  @todo check
  const Vector2D
  localPolarToCartesian(const Vector2D& lpolar) const;

  /// Special method for Disc surface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lcart is local 2D position in carthesian coordinates
  ///
  /// @return value is a local position in polar coordinates
  const Vector2D
  localCartesianToPolar(const Vector2D& lcart) const;

  /// Special method for DiscSurface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lpolar is a local position in polar coordinates
  ///
  /// @return values is local 2D position in carthesian coordinates
  const Vector2D
  localPolarToLocalCartesian(const Vector2D& lpolar) const;

  /// Special method for DiscSurface :  local<->global transformation when
  /// provided cartesian coordinates
  ///
  /// @param lcart is local 2D position in carthesian coordinates
  ///
  /// @return value is a global carthesian 3D position
  const Vector3D
  localCartesianToGlobal(const Vector2D& lcart) const;

  /// Special method for DiscSurface : global<->local from cartesian coordinates
  ///
  /// @param gpos is a global carthesian 3D position
  /// @param tol is the absoltue tolerance parameter
  ///
  /// @return value is a local polar
  const Vector2D
  globalToLocalCartesian(const Vector3D& gpos, double tol = 0.) const;

  /// Path correction due to incident of the track
  ///
  /// @param gpos is the global position as a starting point
  /// @param mom is the global momentum at the starting point
  ///
  /// @return is the correction factor due to incident
  double
  pathCorrection(const Vector3D& gpos,
                 const Vector3D& mom) const final override;

  /// fast straight line intersection schema - standard: provides closest
  /// intersection and (signed) path length
  ///  forceDir is to provide the closest forward solution
  ///
  /// @param gpos is the global position as a starting point
  /// @param dir is the global direction at the starting point
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
  ///  - perpenticular to the normal of the plane
  ///
  /// @return is the intersection object
  virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      dir,
                       bool                 forceDir = false,
                       const BoundaryCheck& bcheck
                       = false) const final override;

  /// Return properly formatted class name for screen output
  virtual std::string
  name() const override
  {
    return "Acts::DiscSurface";
  }

protected:
  std::shared_ptr<const DiscBounds> m_bounds;  ///< bounds (shared)
};

inline DiscSurface*
DiscSurface::clone(const Transform3D* shift) const
{
  if (shift) return new DiscSurface(*this, *shift);
  return new DiscSurface(*this);
}

inline const SurfaceBounds&
DiscSurface::bounds() const
{
  if (m_bounds) return (*(m_bounds.get()));
  return s_noBounds;
}

inline const Vector3D
DiscSurface::normal(const Vector2D&) const
{
  // fast access via tranform matrix (and not rotation())
  auto tMatrix = transform().matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline const Vector3D DiscSurface::binningPosition(BinningValue) const
{
  return center();
}

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

inline double
DiscSurface::pathCorrection(const Vector3D&, const Vector3D& mom) const
{
  /// we can ignore the global position here
  return 1. / std::abs(normal().dot(mom.unit()));
}

inline Intersection
DiscSurface::intersectionEstimate(const Vector3D&      gpos,
                                  const Vector3D&      dir,
                                  bool                 forceDir,
                                  const BoundaryCheck& bcheck) const
{
  double denom = dir.dot(normal());
  if (denom) {
    double   u = (normal().dot((center() - gpos))) / (denom);
    Vector3D intersectPoint(gpos + u * dir);
    // evaluate the intersection in terms of direction
    bool isValid = forceDir ? (u > 0.) : true;
    // evaluate (if necessary in terms of boundaries)
    isValid
        = bcheck ? (isValid && isOnSurface(intersectPoint, bcheck)) : isValid;
    // return the result
    return Intersection(intersectPoint, u, isValid);
  }
  return Intersection(gpos, 0., false);
}

}  // end of namespace

#endif  // ACTS_SURFACES_DISCSURFACE_H
