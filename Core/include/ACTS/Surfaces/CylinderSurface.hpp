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

/// @class CylinderSurface
///
/// Class for a CylinderSurface in the TrackingGeometry.
/// It inherits from Surface.
///
/// The cylinder surface has a special role in the TrackingGeometry,
/// since it builds the surfaces of all TrackingVolumes at container level
/// for a cylindrical tracking geometry.
///
/// @image html CylinderSurface.gif

class CylinderSurface : public Surface
{
public:
  /// Default Constructor - deleted*/
  CylinderSurface() = delete;

  /// Constructor from Transform3D, radius and halflenght
  ///
  /// @param htrans transform to position the surface, can be nullptr
  /// @note if htrans == nullptr, the cylinder is positioned around (0.,0.,0.)
  /// @param radius is the radius of the cylinder
  /// @param hlength is the half lenght of the cylinder in z
  CylinderSurface(std::shared_ptr<Transform3D> htrans,
                  double                       radius,
                  double                       hlength);

  /// Constructor from Transform3D, radius halfphi, and halflenght
  ///
  /// @param htrans transform to position the surface, can be nullptr
  /// @note if htrans == nullptr, the cylinder is positioned around (0.,0.,0.)
  /// @param radius is the radius of the cylinder
  /// @param hphi is the half length in phi of the cylinder
  /// @param hlength is the half lenght of the cylinder in z
  CylinderSurface(std::shared_ptr<Transform3D> htrans,
                  double                       radius,
                  double                       hphi,
                  double                       hlength);

  /// Constructor from Transform3D and CylinderBounds
  ///
  /// @param htrans transform to position the surface, can be nullptr
  /// @note if htrans == nullptr, the cylinder is positioned around (0.,0.,0.)
  /// @param cbounds is a shared pointer to a cylindeer bounds object, must
  /// exist
  CylinderSurface(std::shared_ptr<Transform3D>          htrans,
                  std::shared_ptr<const CylinderBounds> cbounds);

  /// Copy constructor
  ///
  /// @param csf is the source cylinder for the copy
  CylinderSurface(const CylinderSurface& csf);

  /// Copy constructor with shift
  ///
  /// @param csf is the source cylinder for the copy
  /// @param htrans is the additional transform applied after copying the
  /// cylinder
  CylinderSurface(const CylinderSurface& csf, const Transform3D& htrans);

  /// Destructor
  virtual ~CylinderSurface();

  /// Assignment operator
  ///
  /// @param csf is the source cylinder for the copy
  CylinderSurface&
  operator=(const CylinderSurface& csf);

  /// Implicit Constructor - optionally with a shift
  ///
  /// @param shift is an optional transform for a shift applied after coping
  virtual CylinderSurface*
  clone(const Transform3D* shift = nullptr) const override;

  /// The binning position method - is overloaded for r-type binning
  ///
  /// @param bValue is the type of global binning to be done
  ///
  /// @return is the global position to be used for binning
  virtual const Vector3D
  binningPosition(BinningValue bValue) const override;

  /// Return the measurement frame - this is needed for alignment, in particular
  /// The measurement frame of a cylinder is the tangential plane at a given
  /// position
  ///
  /// @param gpos is the position where the measurement frame is defined
  /// @param mom is the momentum vector (ignored)
  ///
  /// @return rotation matrix that defines the measurement frame
  virtual const RotationMatrix3D
  measurementFrame(const Vector3D& gpos, const Vector3D& mom) const override;

  /// Return the surface type
  virtual SurfaceType
  type() const override
  {
    return Surface::Cylinder;
  }

  /// Return method for surface normal information
  /// @note for a Cylinder a local position is always required for the normal
  /// vector
  ///
  /// @param lpos is the local postion for which the normal vector is requested
  ///
  /// @return normal vector at the local position
  virtual const Vector3D
  normal(const Vector2D& lpos) const override;

  /// Return method for surface normal information
  /// @note for a Cylinder a local position is always required for the normal
  /// vector
  ///
  /// @param gpos is the global postion for which the normal vector is requested
  ///
  /// @return normal vector at the global position
  virtual const Vector3D
  normal(const Vector3D& gpos) const override;

  /// Return method for the rotational symmetry axis
  /// @return  the z-Axis of transform
  virtual const Vector3D
  rotSymmetryAxis() const;

  /// This method returns the CylinderBounds by reference
  virtual const CylinderBounds&
  bounds() const override;

  /// Local to global transformation
  ///
  /// @param lpos is the local position to be transformed
  /// @param mom is the global momentum (ignored in this operation)
  /// @param gpos is the global position shich is filled
  virtual void
  localToGlobal(const Vector2D& lpos,
                const Vector3D& mom,
                Vector3D&       gpos) const override;

  /// Global to local transfomration
  ///
  /// @param gpos is the global position to be transformed
  /// @param mom is the global momentum (ignored in this operation)
  /// @param lpos is hte local position to be filled
  ///
  /// @return is a boolean indicating if the transformation succeeded
  virtual bool
  globalToLocal(const Vector3D& gpos,
                const Vector3D& mom,
                Vector2D&       lpos) const override;

  /// Check for position on surface
  ///
  /// @param gpos is the global position to be checked
  /// @param bcheck is the boundary check object
  /// @return is a boolean indicating if the position is on surface
  virtual bool
  isOnSurface(const Vector3D&      gpos,
              const BoundaryCheck& bcheck = true) const override;

  /// Fast straight line intersection schema - provides closest intersection
  ///  and (signed) path length
  ///
  /// @param gpos is the global position as a starting point
  /// @param dir is the global direction at the starting point
  /// @param forceDir is a boolean forcing a solution along direction
  /// @param bchk is the boundary check
  ///
  ///  <b>mathematical motivation:</b>
  ///
  ///  The calculation will be done in the 3-dim frame of the cylinder,
  ///  i.e. the symmetry axis of the cylinder is the z-axis, x- and y-axis are
  /// perpenticular
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
                       const Vector3D&      dir,
                       bool                 forceDir = false,
                       const BoundaryCheck& bchk     = false) const override;

  /// Path correction due to incident of the track
  ///
  /// @param gpos is the global position as a starting point
  /// @param mom is the global momentum at the starting point
  ///
  /// @return is the correction factor due to incident
  virtual double
  pathCorrection(const Vector3D& gpos, const Vector3D& mom) const override;

  /// Return method for properly formatted output string
  virtual std::string
  name() const override
  {
    return "Acts::CylinderSurface";
  }

protected:
  mutable std::shared_ptr<const CylinderBounds> m_bounds;  //!< bounds (shared)
};

inline CylinderSurface*
CylinderSurface::clone(const Transform3D* shift) const
{
  if (shift) return new CylinderSurface(*this, *shift);
  return new CylinderSurface(*this);
}

inline const Vector3D
CylinderSurface::normal(const Vector2D& lpos) const
{
  double   phi = lpos[Acts::eLOC_RPHI] / m_bounds->r();
  Vector3D localNormal(cos(phi), sin(phi), 0.);
  return Vector3D(transform().rotation() * localNormal);
}

inline const Vector3D
CylinderSurface::normal(const Vector3D& gpos) const
{
  // get it into the cylinder frame if needed
  Vector3D pos3D = gpos;
  if (m_transform || m_associatedDetElement) {
    pos3D = transform().inverse() * gpos;
  }
  // set the z coordinate to 0
  pos3D.z() = 0.;
  return pos3D.unit();
}

inline double
CylinderSurface::pathCorrection(const Vector3D& gpos, const Vector3D& mom) const
{
  Vector3D normalT  = normal(gpos);
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
