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
  /// @param htrans transform to position the surface, can be nullptr
  /// @note if htrans == nullptr, the cylinder is positioned around (0.,0.,0.)
  /// @param radius is the radius of the cylinder
  /// @param hlength is the half lenght of the cylinder in z 
  CylinderSurface(std::shared_ptr<Transform3D> htrans,
                  double                       radius,
                  double                       hlength);

  /// Constructor from Transform3D, radius halfphi, and halflenght
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
  /// @param htrans transform to position the surface, can be nullptr
  /// @note if htrans == nullptr, the cylinder is positioned around (0.,0.,0.)
  /// @param cbounds is a shared pointer to a cylindeer bounds object                
  CylinderSurface(std::shared_ptr<Transform3D>          htrans,
                  std::shared_ptr<const CylinderBounds> cbounds);

  /// Copy constructor 
  /// @param csf is the source cylinder for the copy                  
  CylinderSurface(const CylinderSurface& csf);

  /// Copy constructor with shift
  /// @param csf is the source cylinder for the copy
  /// @param htrans is the additional transform applied after copying the cylinder                  
  CylinderSurface(const CylinderSurface& csf, const Transform3D& htrans);

  /// Destructor
  virtual ~CylinderSurface();

  /// Assignment operator
  CylinderSurface&
  operator=(const CylinderSurface& csf);

  /// Implicit Constructor - optionally with a shift 
  /// @param shift is an optional transform for a shift applied after coping
  virtual CylinderSurface*
  clone(const Transform3D* shift = nullptr) const override;

  /// The binning position method - is overloaded for r-type binning 
  /// @param bValue is the type of global binning to be done
  virtual Vector3D
  binningPosition(BinningValue bValue) const override;

  /// Return the measurement frame - this is needed for alignment, in particular
  /// The measurement frame of a cylinder is the tangential plane at a given position
  /// @param gpos is the position where the measurement frame is defined
  /// @param mom is the momentum vector (ignored)
  virtual const RotationMatrix3D
  measurementFrame(const Vector3D& gpos,
                   const Vector3D& mom) const override;

  /// Return the surface type 
  virtual SurfaceType
  type() const override
  {
    return Surface::Cylinder;
  }

  /// Return method for surface normal information
  /// @note for a Cylinder a local position is always required for the normal vector
  /// @param lpos is the local postion for which the normal vector is requested
  virtual const Vector3D
  normal(const Vector2D& lpos) const override;

  ///  Return method for the rotational symmetry axis - the z-Axis of the HepTransform
  virtual const Vector3D
  rotSymmetryAxis() const;

  /// This method returns the CylinderBounds by reference
  virtual const CylinderBounds&
  bounds() const override;

  /// @copydoc Surface::localToGlobal
  /// @note momentum is ignored int he calculation 
  virtual void
  localToGlobal(const Vector2D& lpos,
                const Vector3D& mom,
                Vector3D&       gpos) const override;

  /// @copydoc Surface::globalToLocal
  /// @note momentum is ignored int he calculation 
  virtual bool
  globalToLocal(const Vector3D& gpos,
                const Vector3D& mom,
                Vector2D&       lpos) const override;

  /// @copydoc Surface::isOnSurface
  virtual bool
  isOnSurface(const Vector3D&      gpos,
              const BoundaryCheck& bchk = true) const override;

  ///  fast straight line intersection schema - provides closest intersection and
  /// (signed) path length
  /// 
  /// @copydoc Surface::intersectionEstimate
  /// 
  ///  <b>mathematical motivation:</b>
  /// 
  ///  The calculation will be done in the 3-dim frame of the cylinder,
  ///  i.e. the symmetry axis of the cylinder is the z-axis, x- and y-axis are
  /// perpenticular
  ///  to the the z-axis. In this frame the cylinder is centered around the
  /// origin.
  ///  Therefore the two points describing the line have to be first recalculated
  /// into the new frame.
  ///  Suppose, this is done, the intersection is straight forward:<br>
  ///  may @f$p_{1}=(p_{1x}, p_{1y}, p_{1z}), p_{2}=(p_{2x}, p_{2y}, p_{2z})
  /// @f$the two points describing the 3D-line,
  ///  then the line in the \f$x-y@f$plane can be written as
  ///  @f$y=kx+d\f$, where @f$k =\frac{p_{2y}-p_{1y}}{p_{2x}-p_{1x}}@f$such as
  /// @f$d=\frac{p_{2x}p_{1y}-p_{1x}p_{2y}}{p_{2x}-p_{1x}},\f$<br>
  ///  and intersects with the corresponding circle @f$x^{2}+y^{2} = R^{2}.
  /// @f$<br>
  ///  The solutions can then be found by a simple quadratic equation and
  /// reinsertion into the line equation.
  /// 
  virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      dir,
                       bool                 forceDir = false,
                       const BoundaryCheck& bchk     = false) const override;

  /// Path correction due to incident
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
CylinderSurface::normal(const Vector2D& lp) const
{
  double   phi = lp[Acts::eLOC_RPHI] / (bounds().r());
  Vector3D localNormal(cos(phi), sin(phi), 0.);
  return std::move(Vector3D(transform().rotation() * localNormal));
}

inline double
CylinderSurface::pathCorrection(const Vector3D& gpos, const Vector3D& mom) const
{
  // the global normal vector is pos-center.unit() - at the z of the position
  if (m_transform || m_associatedDetElement){
    Transform3D iTransform = transform().inverse();
    return pathCorrection(iTransform*gpos,iTransform.linear()*gpos);
  }
  // 
  Vector3D pcT(gpos.x() - center().x(), gpos.y() - center().y(), 0.);
  Vector3D normalT(pcT.unit());  // transverse normal
  double   cosAlpha = normalT.dot(mom.unit());
  return fabs(1./cosAlpha);
}

inline const CylinderBounds&
CylinderSurface::bounds() const
{
  return (*m_bounds.get());
}

}  // end of namespace

#endif  // ACTS_SURFACES_CYLINDERSURFACE_H
