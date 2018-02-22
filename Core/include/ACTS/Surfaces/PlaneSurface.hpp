// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlaneSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_PLANESURFACE_H
#define ACTS_SURFACES_PLANESURFACE_H 1

#include <limits>
#include "ACTS/Surfaces/InfiniteBounds.hpp"
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Identifier.hpp"

namespace Acts {

class DetectorElementBase;

/// @class PlaneSurface
///
/// Class for a planaer in the TrackingGeometry.
///
/// The PlaneSurface extends the Surface class with the possibility to
/// convert local to global positions (vice versa).
///
/// @image html PlaneSurface.png
///
class PlaneSurface : public Surface
{
public:
  /// Default Constructor - needed for persistency
  PlaneSurface();

  /// Copy Constructor
  ///
  /// @param psf is the source surface for the copy
  PlaneSurface(const PlaneSurface& other);

  /// Copy Constructor with shift
  ///
  /// @param psf is the source surface for the copy
  /// @param htrans is the transformation that positions the surface in space
  PlaneSurface(const PlaneSurface& other, const Transform3D& htrans);

  /// Dedicated Constructor with normal vector
  /// This is for curvilinear surfaces which are by definition boundless
  ///
  /// @param center is the center position of the surface
  /// @param normal is thenormal vector of the plane surface
  PlaneSurface(const Vector3D& center, const Vector3D& normal);

  /// Constructor from DetectorElementBase
  ///
  /// @param pbounds are the provided planar bounds (shared)
  /// @param detelement is the linked detector element to this surface
  /// @param identifier is the identifier of associated to this surfacex
  PlaneSurface(std::shared_ptr<const PlanarBounds> pbounds,
               const DetectorElementBase&          detelement,
               const Identifier&                   identifier = Identifier());

  /// Constructor for Planes with (optional) shared bounds object
  ///
  /// @param htrans transform in 3D that positions this surface
  /// @param pbounds bounds object to describe the actual surface area
  PlaneSurface(std::shared_ptr<const Transform3D>  htrans,
               std::shared_ptr<const PlanarBounds> pbounds = nullptr);

  virtual ~PlaneSurface();

  /// Assignment operator
  ///
  /// @param psf source PlaneSurface for assignment
  PlaneSurface&
  operator=(const PlaneSurface& other);

  /// Virtual constructor with optional shift
  /// ownership of the shift transform is not given !!
  ///
  /// @param shift is a potential shift after cloning
  virtual PlaneSurface*
  clone(const Transform3D* shift = nullptr) const override;

  /// Normal vector return
  ///
  /// @param lpos is the local position is ignored
  /// return a Vector3D by value
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

  /// Return the surface type
  virtual SurfaceType
  type() const override;

  /// Return method for bounds object of this surfrace
  virtual const SurfaceBounds&
  bounds() const override;

  /// Geometrical on surface test
  /// This method returns true if the GlobalPosition is on the Surface for both,
  /// within or without check of whether the local position is inside boundaries
  /// or not
  ///
  /// @param gpos global position to be checked
  /// @param bcheck gboundary check directive
  ///
  /// @return is a boolean indicator if the position is on surface
  virtual bool
  isOnSurface(const Vector3D&      gpos,
              const BoundaryCheck& bcheck = true) const override;

  /// Local to global transformation
  /// For planar surfaces the momentum is ignroed in the local to global
  /// transformation
  ///
  /// @param lpos local 2D posittion in specialized surface frame
  /// @param mom global 3D momentum representation (optionally ignored)
  /// @param gpos global 3D position to be filled (given by reference for method
  /// symmetry)
  virtual void
  localToGlobal(const Vector2D& lpos,
                const Vector3D& mom,
                Vector3D&       gpos) const override;

  /// Global to local transformation
  /// For planar surfaces the momentum is ignroed in the global to local
  /// transformation
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
                Vector2D&       lpos) const override;

  /// Method that calculates the correction due to incident angle
  ///
  /// @param gpos global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param mom global 3D momentum representation (optionally ignored)
  ///
  /// @note this is the final implementation of the pathCorrection function
  ///
  /// @return a double representing the scaling factor
  double
  pathCorrection(const Vector3D& gpos, const Vector3D& mom) const final;

  ///  fast straight line intersection schema - standard: provides closest
  /// intersection and (signed) path length
  ///  forceDir is to provide the closest forward solution
  ///
  ///  @param gpos is the start position of the intersection attempt
  ///  @param gdir is the direction of the interesection attempt,
  ///        @note has to be normalized
  ///  @param forceDir is the directive whether to force only foward solution
  ///  (w.r.t dir)
  ///  @param bcheck is the boundary check directive
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
  virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      gdir,
                       bool                 forceDir = true,
                       const BoundaryCheck& bcheck   = true) const override;

  /// Return properly formatted class name for screen output
  virtual std::string
  name() const override;

protected:
  /// the bounds of this surface
  std::shared_ptr<const PlanarBounds> m_bounds;
};

inline Intersection
PlaneSurface::intersectionEstimate(const Vector3D&      gpos,
                                   const Vector3D&      gdir,
                                   bool                 forceDir,
                                   const BoundaryCheck& bcheck) const
{
  double denom = gdir.dot(normal());
  if (denom) {
    double   u = (normal().dot((center() - gpos))) / (denom);
    Vector3D intersectPoint(gpos + u * gdir);
    // evaluate the intersection in terms of direction
    bool isValid = forceDir ? (u > 0.) : true;
    // evaluate (if necessary in terms of boundaries)
    isValid
        = bcheck ? (isValid && isOnSurface(intersectPoint, bcheck)) : isValid;
    // return the result
    return Intersection(intersectPoint, u, isValid);
  }
  return Intersection(gpos, std::numeric_limits<double>::max(), false);
}

}  // end of namespace

#endif  // ACTS_SURFACES_PLANESURFACE_H
