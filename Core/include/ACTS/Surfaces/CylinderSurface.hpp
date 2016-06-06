// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_CYLINDERSURFACE_H
#define ACTS_SURFACES_CYLINDERSURFACE_H 1

#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

/**
 @class CylinderSurface

 Class for a CylinderSurface in the TrackingGeometry.
 It inherits from Surface.

 The cylinder surface has a special role in the TrackingGeometry,
 since it builds the surfaces of all TrackingVolumes at container level
 for a cylindrical tracking geometry.

 @image html CylinderSurface.gif

 */

class CylinderSurface : public Surface
{
public:
  /** Default Constructor*/
  CylinderSurface();

  /** Constructor from Transform3D, radius and halflenght*/
  CylinderSurface(std::shared_ptr<Transform3D> htrans,
                  double                       radius,
                  double                       hlength);

  /** Constructor from Transform3D, radius halfphi, and halflenght*/
  CylinderSurface(std::shared_ptr<Transform3D> htrans,
                  double                       radius,
                  double                       hphi,
                  double                       hlength);

  /** Constructor from Transform3D and CylinderBounds
    - ownership of the bounds is passed */
  CylinderSurface(std::shared_ptr<Transform3D> htrans,
                  const CylinderBounds*        cbounds);

  /** Constructor from Transform3D and CylinderBounds
    - ownership of the bounds is passed */
  CylinderSurface(std::shared_ptr<Transform3D>          htrans,
                  std::shared_ptr<const CylinderBounds> cbounds);

  /** Constructor from Transform3D from unique_ptr.
     - bounds is not set */
  CylinderSurface(std::unique_ptr<Transform3D> htrans);

  /** Constructor from radius and halflenght - speed optimized for concentric
   * volumes */
  CylinderSurface(double radius, double hlength);

  /** Constructor from radius halfphi, and halflenght - speed optimized fron
   * concentric volumes */
  CylinderSurface(double radius, double hphi, double hlength);

  /** Copy constructor */
  CylinderSurface(const CylinderSurface& csf);

  /** Copy constructor with shift */
  CylinderSurface(const CylinderSurface& csf, const Transform3D& transf);

  /** Destructor*/
  virtual ~CylinderSurface();

  /** Assignment operator*/
  CylinderSurface&
  operator=(const CylinderSurface& csf);

  /** Equality operator*/
  virtual bool
  operator==(const Surface& sf) const override;

  /** Implicit Constructor - optionally with a shift */
  virtual CylinderSurface*
  clone(const Transform3D* shift = nullptr) const override;

  /** The binning position method - is overloaded for r-type binning */
  virtual Vector3D
  binningPosition(BinningValue bValue) const override;

  /** Return the measurement frame - this is needed for alignment, in particular
     for StraightLine and Perigee Surface
      - the default implementation is the the RotationMatrix3D of the transform
     */
  virtual const RotationMatrix3D
  measurementFrame(const Vector3D& glopos,
                   const Vector3D& glomom) const override;

  /** Return the surface type */
  virtual SurfaceType
  type() const override
  {
    return Surface::Cylinder;
  }

  /** Return method for surface normal information
     at a given local point, overwrites the normal() from base class.*/
  virtual const Vector3D&
  normal() const override;

  /** Return method for surface normal information
     at a given local point, overwrites the normal() from base class.*/
  virtual const Vector3D
  normal(const Vector2D& locpo) const override;

  /** Return method for surface normal information
     at a given global point, overwrites the normal() from base class.*/
  virtual const Vector3D
  normal(const Vector3D& global) const override;

  /** Return method for the rotational symmetry axis - the z-Axis of the
   * HepTransform */
  virtual const Vector3D&
  rotSymmetryAxis() const;

  /** This method returns the CylinderBounds by reference
   (NoBounds is not possible for cylinder)*/
  virtual const CylinderBounds&
  bounds() const override;

  /** Specialized for CylinderSurface : LocalToGlobal method without dynamic
   * memory allocation */
  virtual void
  localToGlobal(const Vector2D& locp,
                const Vector3D& mom,
                Vector3D&       glob) const override;

  /** Specialized for CylinderSurface : GlobalToLocal method without dynamic
   * memory allocation - boolean checks if on surface */
  virtual bool
  globalToLocal(const Vector3D& glob,
                const Vector3D& mom,
                Vector2D&       loc) const override;

  /** This method returns true if the GlobalPosition is on the Surface for both,
    within
    or without check of whether the local position is inside boundaries or not
    */
  virtual bool
  isOnSurface(const Vector3D&      glopo,
              const BoundaryCheck& bchk = true) const override;

  /** fast straight line intersection schema - provides closest intersection and
     (signed) path length

      <b>mathematical motivation:</b>

      The calculation will be done in the 3-dim frame of the cylinder,
      i.e. the symmetry axis of the cylinder is the z-axis, x- and y-axis are
     perpenticular
      to the the z-axis. In this frame the cylinder is centered around the
     origin.
      Therefore the two points describing the line have to be first recalculated
     into the new frame.
      Suppose, this is done, the intersection is straight forward:<br>
      may @f$p_{1}=(p_{1x}, p_{1y}, p_{1z}), p_{2}=(p_{2x}, p_{2y}, p_{2z})
     @f$the two points describing the 3D-line,
      then the line in the \f$x-y@f$plane can be written as
      @f$y=kx+d\f$, where @f$k =\frac{p_{2y}-p_{1y}}{p_{2x}-p_{1x}}@f$such as
     @f$d=\frac{p_{2x}p_{1y}-p_{1x}p_{2y}}{p_{2x}-p_{1x}},\f$<br>
      and intersects with the corresponding circle @f$x^{2}+y^{2} = R^{2}.
     @f$<br>
      The solutions can then be found by a simple quadratic equation and
     reinsertion into the line equation.
  */
  virtual Intersection
  intersectionEstimate(const Vector3D&      pos,
                       const Vector3D&      dir,
                       bool                 forceDir = false,
                       const BoundaryCheck& bchk     = false) const override;

  /** the pathCorrection for derived classes with thickness */
  virtual double
  pathCorrection(const Vector3D& pos, const Vector3D& mom) const override;

  /** Return properly formatted class name for screen output */
  virtual std::string
  name() const override
  {
    return "Acts::CylinderSurface";
  }

protected:                                                 //!< data members
  mutable std::shared_ptr<const CylinderBounds> m_bounds;  //!< bounds (shared)
  mutable Vector3D* m_rotSymmetryAxis;  //!< The rotational symmetry axis
};

inline CylinderSurface*
CylinderSurface::clone(const Transform3D* shift) const
{
  if (shift) return new CylinderSurface(*this, *shift);
  return new CylinderSurface(*this);
}

inline const Vector3D&
CylinderSurface::normal() const
{
  return Surface::normal();
}

inline const Vector3D
CylinderSurface::normal(const Vector2D& lp) const
{
  double   phi = lp[Acts::eLOC_RPHI] / (bounds().r());
  Vector3D localNormal(cos(phi), sin(phi), 0.);
  return Vector3D(transform().rotation() * localNormal);
}

inline const Vector3D
CylinderSurface::normal(const Vector3D& global) const
{
  // transform the global position into the the cylinder frame
  const Vector3D& local3D = m_transform
      ? Vector3D(transform().inverse().linear() * global)
      : global;
  double   phi = local3D.phi();
  Vector3D localNormal(cos(phi), sin(phi), 0.);
  return m_transform ? Vector3D(transform().rotation() * localNormal)
                     : localNormal;
}

inline double
CylinderSurface::pathCorrection(const Vector3D& pos, const Vector3D& mom) const
{
  // the global normal vector is pos-center.unit() - at the z of the position
  // @TODO make safe for tilt
  Vector3D pcT(pos.x() - center().x(), pos.y() - center().y(), 0.);
  Vector3D normalT(pcT.unit());  // transverse normal
  double   cosAlpha = normalT.dot(mom.unit());
  return fabs(1. / cosAlpha);
}

inline const CylinderBounds&
CylinderSurface::bounds() const
{
  return (*m_bounds.get());
}

}  // end of namespace

#endif  // ACTS_SURFACES_CYLINDERSURFACE_H
