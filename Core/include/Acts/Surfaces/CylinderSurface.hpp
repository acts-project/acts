// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderSurface.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <cmath>

#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PolyhedronRepresentation.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"

namespace Acts {

/// @class CylinderSurface
///
/// Class for a CylinderSurface in the TrackingGeometry.
/// It inherits from Surface.
///
/// The cylinder surface has a special role in the TrackingGeometry,
/// since it builds the surfaces of all TrackingVolumes at container level
/// for a cylindrical tracking geometry.
///
/// @image html CylinderSurface.png

class CylinderSurface : public Surface
{
public:
  /// Deleted default constructor
  CylinderSurface() = delete;

  /// Constructor from Transform3D, radius and halflenght
  ///
  /// @param htrans transform to position the surface, can be nullptr
  /// @note if htrans == nullptr, the cylinder is positioned around (0.,0.,0.)
  /// @param radius is the radius of the cylinder
  /// @param hlength is the half lenght of the cylinder in z
  CylinderSurface(std::shared_ptr<const Transform3D> htrans,
                  double                             radius,
                  double                             hlength);

  /// Constructor from Transform3D, radius halfphi, and halflenght
  ///
  /// @param htrans transform to position the surface, can be nullptr
  /// @note if htrans == nullptr, the cylinder is positioned around (0.,0.,0.)
  /// @param radius is the radius of the cylinder
  /// @param hphi is the half length in phi of the cylinder
  /// @param hlength is the half lenght of the cylinder in z
  CylinderSurface(std::shared_ptr<const Transform3D> htrans,
                  double                             radius,
                  double                             hphi,
                  double                             hlength);

  /// Constructor from DetectorElementBase
  ///
  /// @param cbounds are the provided cylinder bounds (shared)
  /// @param detelement is the linked detector element to this surface
  /// @param identifier is the identifier of associated to this surfacex
  CylinderSurface(std::shared_ptr<const CylinderBounds> cbounds,
                  const DetectorElementBase&            detelement,
                  const Identifier& identifier = Identifier());

  /// Constructor from Transform3D and CylinderBounds
  ///
  /// @param htrans transform to position the surface, can be nullptr
  /// @note if htrans == nullptr, the cylinder is positioned around (0.,0.,0.)
  /// @param cbounds is a shared pointer to a cylindeer bounds object,
  /// it must exist (assert test)
  CylinderSurface(std::shared_ptr<const Transform3D>    htrans,
                  std::shared_ptr<const CylinderBounds> cbounds);

  /// Copy constructor
  ///
  /// @param other is the source cylinder for the copy
  CylinderSurface(const CylinderSurface& other);

  /// Copy constructor with shift
  ///
  /// @param other is the source cylinder for the copy
  /// @param htrans is the additional transform applied after copying the
  /// cylinder
  CylinderSurface(const CylinderSurface& other, const Transform3D& htrans);

  /// Constructor which accepts @c variant_data
  ///
  /// @param data the @c variant_data to build from
  CylinderSurface(const variant_data& data);

  /// Destructor
  virtual ~CylinderSurface();

  /// Assignment operator
  ///
  /// @param other is the source cylinder for the copy
  CylinderSurface&
  operator=(const CylinderSurface& other);

  /// Implicit Constructor - optionally with a shift
  ///
  /// @param shift is an optional transform for a shift applied after coping
  virtual CylinderSurface*
  clone(const Transform3D* shift = nullptr) const final override;

  /// The binning position method - is overloaded for r-type binning
  ///
  /// @param bValue is the type of global binning to be done
  ///
  /// @return is the global position to be used for binning
  virtual const Vector3D
  binningPosition(BinningValue bValue) const final override;

  /// Return the measurement frame - this is needed for alignment, in particular
  /// The measurement frame of a cylinder is the tangential plane at a given
  /// position
  ///
  /// @param gpos is the position where the measurement frame is defined
  /// @param mom is the momentum vector (ignored)
  /// @return rotation matrix that defines the measurement frame
  virtual const RotationMatrix3D
  referenceFrame(const Vector3D& gpos,
                 const Vector3D& mom) const final override;

  /// Return the surface type
  virtual SurfaceType
  type() const override;

  /// Return method for surface normal information
  /// @note for a Cylinder a local position is always required for the normal
  /// vector
  ///
  /// @param lpos is the local postion for which the normal vector is requested
  /// @return normal vector at the local position
  virtual const Vector3D
  normal(const Vector2D& lpos) const final override;

  /// Return method for surface normal information
  /// @note for a Cylinder a local position is always required for the normal
  /// vector
  ///
  /// @param gpos is the global postion for which the normal vector is requested
  /// @return normal vector at the global position
  virtual const Vector3D
  normal(const Vector3D& gpos) const final override;

  /// Return method for the rotational symmetry axis
  ///
  /// @return  the z-Axis of transform
  virtual const Vector3D
  rotSymmetryAxis() const;

  /// This method returns the CylinderBounds by reference
  virtual const CylinderBounds&
  bounds() const final override;

  /// Local to global transformation
  ///
  /// @param lpos is the local position to be transformed
  /// @param mom is the global momentum (ignored in this operation)
  /// @param gpos is the global position shich is filled
  virtual void
  localToGlobal(const Vector2D& lpos,
                const Vector3D& mom,
                Vector3D&       gpos) const final override;

  /// Global to local transfomration
  ///
  /// @param gpos is the global position to be transformed
  /// @param mom is the global momentum (ignored in this operation)
  /// @param lpos is hte local position to be filled
  /// @return is a boolean indicating if the transformation succeeded
  virtual bool
  globalToLocal(const Vector3D& gpos,
                const Vector3D& mom,
                Vector2D&       lpos) const final override;

  /// Check for position on surface
  ///
  /// @param gpos is the global position to be checked
  /// @param bcheck is the boundary check object
  /// @return is a boolean indicating if the position is on surface
  virtual bool
  isOnSurface(const Vector3D&      gpos,
              const BoundaryCheck& bcheck = true) const final override;

  /// Fast straight line intersection schema - provides closest intersection
  ///  and (signed) path length
  ///
  /// @param gpos is the global position as a starting point
  /// @param gdir is the global direction at the starting point,
  ///        @note has to be normalized
  /// @param forceDir is a boolean forcing a solution along direction
  /// @param bcheck is the boundary check
  ///
  ///  <b>mathematical motivation:</b>
  ///
  ///  The calculation will be done in the 3-dim frame of the cylinder,
  ///  i.e. the symmetry axis of the cylinder is the z-axis, x- and y-axis are
  /// perpendicular
  ///  to the the z-axis. In this frame the cylinder is centered around the
  /// origin.
  ///  Therefore the two points describing the line have to be first
  ///  recalculated
  /// into the new frame.
  ///  Suppose, this is done, the intersection is straight forward:
  ///  @f$p_{1}=(p_{1x}, p_{1y}, p_{1z}), p_{2}=(p_{2x}, p_{2y}, p_{2z}) @f$
  ///  are the two points describing the 3D-line,
  ///  then the line in the \f$x-y@f$ plane can be written as
  ///  @f$y=kx+d\f$, where @f$k =\frac{p_{2y}-p_{1y}}{p_{2x}-p_{1x}}@f$such as
  /// @f$d=\frac{p_{2x}p_{1y}-p_{1x}p_{2y}}{p_{2x}-p_{1x}},\f$<br>
  ///  and intersects with the corresponding circle @f$x^{2}+y^{2} = R^{2}.
  /// @f$<br>
  ///  The solutions can then be found by a simple quadratic equation and
  /// reinsertion into the line equation.
  ///
  /// @return is the intersection object
  virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      gdir,
                       bool                 forceDir = false,
                       const BoundaryCheck& bcheck
                       = false) const final override;

  /// Path correction due to incident of the track
  ///
  /// @param gpos is the global position as a starting point
  /// @param mom is the global momentum at the starting point
  /// @return is the correction factor due to incident
  virtual double
  pathCorrection(const Vector3D& gpos,
                 const Vector3D& mom) const final override;

  /// Return method for properly formatted output string
  virtual std::string
  name() const override;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  virtual variant_data
  toVariantData() const override;

  /// Return a PolyhedronRepresentation for this object
  /// @param l0div Number of divisions along l0 (phi)
  /// @param l1div Number of divisions along l1 (z)
  virtual PolyhedronRepresentation
  polyhedronRepresentation(size_t l0div = 10, size_t l1div = 1) const;

protected:
  std::shared_ptr<const CylinderBounds> m_bounds;  //!< bounds (shared)
};

inline const Vector3D
CylinderSurface::rotSymmetryAxis() const
{
  // fast access via tranform matrix (and not rotation())
  return transform().matrix().block<3, 1>(0, 2);
}

inline bool
CylinderSurface::isOnSurface(const Vector3D&      gpos,
                             const BoundaryCheck& bcheck) const
{
  Vector3D loc3Dframe
      = Surface::m_transform ? (transform().inverse()) * gpos : gpos;
  return (bcheck ? bounds().inside3D(loc3Dframe, bcheck) : true);
}

inline Intersection
CylinderSurface::intersectionEstimate(const Vector3D&      gpos,
                                      const Vector3D&      gdir,
                                      bool                 forceDir,
                                      const BoundaryCheck& bcheck) const
{
  bool needsTransform
      = (Surface::m_transform || m_associatedDetElement) ? true : false;
  // create the hep points
  Vector3D point1    = gpos;
  Vector3D direction = gdir;
  if (needsTransform) {
    Transform3D invTrans = transform().inverse();
    point1               = invTrans * gpos;
    direction            = invTrans.linear() * gdir;
  }
  Vector3D point2 = point1 + direction;
  // the bounds radius
  double R  = bounds().r();
  double t1 = 0.;
  double t2 = 0.;
  if (direction.x()) {
    // get line and circle constants
    double k = (direction.y()) / (direction.x());
    double d = (point2.x() * point1.y() - point1.x() * point2.y())
        / (point2.x() - point1.x());
    // and solve the qaudratic equation
    detail::RealQuadraticEquation pquad(1 + k * k, 2 * k * d, d * d - R * R);
    if (pquad.solutions == 2) {
      // the solutions in the 3D frame of the cylinder
      t1 = (pquad.first - point1.x()) / direction.x();
      t2 = (pquad.second - point1.x()) / direction.x();
    } else {  // bail out if no solution exists
      return Intersection(gpos, 0., false);
    }
  } else {
    // bail out if no solution exists
    if (!direction.y()) return Intersection(gpos, 0., false);
    // x value ise th one of point1
    // x^2 + y^2 = R^2
    // y = sqrt(R^2-x^2)
    double x     = point1.x();
    double r2mx2 = R * R - x * x;
    // bail out if no solution
    if (r2mx2 < 0.) return Intersection(gpos, 0., false);
    double y = sqrt(r2mx2);
    // assign parameters and solutions
    t1 = (y - point1.y()) / direction.y();
    t2 = (-y - point1.y()) / direction.y();
  }
  Vector3D sol1raw(point1 + t1 * direction);
  Vector3D sol2raw(point1 + t2 * direction);
  // now reorder and return
  Vector3D solution(0, 0, 0);
  double   path = 0.;

  // first check the validity of the direction
  bool isValid = true;

  // both solutions are of same sign, take the smaller, but flag as false if not
  // forward
  if (t1 * t2 > 0 || !forceDir) {
    // asign validity
    isValid = forceDir ? (t1 > 0.) : true;
    // assign the right solution
    if (t1 * t1 < t2 * t2) {
      solution = sol1raw;
      path     = t1;
    } else {
      solution = sol2raw;
      path     = t2;
    }
  } else {
    if (t1 > 0.) {
      solution = sol1raw;
      path     = t1;
    } else {
      solution = sol2raw;
      path     = t2;
    }
  }
  // the solution is still in the local 3D frame, direct check
  isValid = bcheck ? (isValid && bounds().inside3D(solution, bcheck)) : isValid;

  // now return
  return needsTransform ? Intersection((transform() * solution), path, isValid)
                        : Intersection(solution, path, isValid);
}
}  // end of namespace
