// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// StraightLineSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACESSTRAIGHTLINESURFACE_H
#define ACTS_SURFACESSTRAIGHTLINESURFACE_H

#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Surfaces/BoundlessT.hpp"
#include "ACTS/Surfaces/LineBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Identifier.hpp"

namespace Acts {

class DetectorElementBase;

///  @class StraightLineSurface
/// 
///  Class for a StraightLineSurface in the TrackingGeometry
///  to describe dirft tube and straw like detectors.
///  It inherits from Surface.

class StraightLineSurface : public Surface
{
public:
  /// Default Constructor - deleted
  StraightLineSurface() = delete;

  /// Constructor from Transform3D (boundless surface)
  /// @param htrans is the transform that positions the surface in the global frame
  StraightLineSurface(std::shared_ptr<Transform3D> htrans);

  /// Constructor from Transform3D and bounds
  /// @param htrans is the transform that positions the surface in the global frame
  /// @param radius is the straw radius
  /// @param halex is the half length in z
  StraightLineSurface(std::shared_ptr<Transform3D> htrans,
                      double                       radius,
                      double                       halez);

  /// Constructor from Transform3D and a shared bounds object
  /// @param htrans is the transform that positions the surface in the global frame
  /// @param lbounds are teh bounds describing the straw dimensions, can be optionally nullptr                 
  StraightLineSurface(std::shared_ptr<Transform3D>      htrans,
                      std::shared_ptr<const LineBounds> lbounds = nullptr);

  /// Constructor from DetectorElementBase and Element identifier
  /// @param lbounds are teh bounds describing the straw dimensions, they must not be nullptr                    
  /// @param detelement for which this surface is (at least) one representation
  /// @param identifier                                          
  StraightLineSurface(std::shared_ptr<const LineBounds> lbounds,
                      const DetectorElementBase& detelement,
                      const Identifier&          identifier = Identifier());

  /// Copy constructor
  /// @param slsf is teh source surface for copying                      
  StraightLineSurface(const StraightLineSurface& slsf);

  /// Copy constructor with shift
  /// @param slsf is the source surface dor copying
  /// @param transf is the additional transform applied after copying
  StraightLineSurface(const StraightLineSurface& slsf,
                      const Transform3D&         transf);

  /// Destructor
  virtual ~StraightLineSurface();

  /// Assignment operator
  StraightLineSurface&
  operator=(const StraightLineSurface& slsf);

  /// Implicit constructor - shift can be provided */
  virtual StraightLineSurface*
  clone(const Transform3D* shift = nullptr) const override;

  /// Return the measurement frame - this is needed for alignment, in particular
  /// for StraightLine and Perigee Surface
  ///  - the default implementation is the the RotationMatrix3D of the transform
  virtual const RotationMatrix3D
  measurementFrame(const Vector3D& gpos,
                   const Vector3D& mom) const override;

  /// Return the surface type 
  virtual SurfaceType
  type() const override
  {
    return Surface::Line;
  }

  /// @copydoc Surface::localToGlobal
  virtual void
  localToGlobal(const Vector2D& locp,
                const Vector3D& mom,
                Vector3D&       glob) const override;

  /// Specified for StraightLineSurface: GlobalToLocal method without dynamic
  /// memory allocation
  /// This method is the true global->local transformation.<br>
  /// makes use of globalToLocal and indicates the sign of the Acts::eLOC_R by the
  /// given momentum
  /// 
  /// The calculation of the sign of the radius (or \f$ d_0 \f$) can be done as
  /// follows:<br>
  /// May \f$ \vec d = \vec m - \vec c \f$ denote the difference between the
  /// center of the line and
  /// the global position of the measurement/predicted state, then \f$ \vec d \f$
  /// lies within the so
  /// called measurement plane.
  /// The measurement plane is determined by the two orthogonal vectors \f$
  /// \vec{measY}= \vec{Acts::eLOC_Z} \f$
  /// and \f$ \vec{measX} = \vec{measY} \times \frac{\vec{p}}{|\vec{p}|} \f$.<br>
  /// 
  /// The sign of the radius (\f$ d_{0} \f$ ) is then defined by the projection of
  /// \f$ \vec{d} \f$
  /// onto \f$ \vec{measX} \f$:<br>
  /// \f$ sign = -sign(\vec{d} \cdot \vec{measX}) \f$
  /// 
  /// \image html SignOfDriftCircleD0.gif
  virtual bool
  globalToLocal(const Vector3D& glob,
                const Vector3D& mom,
                Vector2D&       loc) const override;

  /// Special method for StraightLineSurface 
  /// provides the Line direction from cache: speedup
  const Vector3D
  lineDirection() const;

  /// fast straight line intersection schema - standard: provides closest
  ///  intersection and (signed) path length
  ///   forceDir is to provide the closest forward solution
  /// 
  ///   b>mathematical motivation:</b>
  ///   Given two lines in parameteric form:<br>
  ///   - @f$ \vec l_{a}(\lambda) = \vec m_a + \lambda \cdot \vec e_{a} @f$ <br>
  ///   - @f$ \vec l_{b}(\mu) = \vec m_b + \mu \cdot \vec e_{b} @f$ <br>
  ///   the vector between any two points on the two lines is given by:
  ///   - @f$ \vec s(\lambda, \mu) = \vec l_{b} - l_{a} = \vec m_{ab} + \mu \cdot
  ///  \vec e_{b} - \lambda \cdot \vec e_{a} @f$, <br>
  ///   when @f$ \vec m_{ab} = \vec m_{b} - \vec m_{a} @f$.<br>
  ///   @f$ \vec s(\lambda_0, \mu_0) @f$  denotes the vector between the two
  ///  closest points <br>
  ///   @f$ \vec l_{a,0} = l_{a}(\lambda_0) @f$ and @f$ \vec l_{b,0} =
  ///  l_{b}(\mu_0) @f$ <br>
  ///   and is perpenticular to both, @f$ \vec e_{a} @f$ and @f$ \vec e_{b} @f$.
  /// 
  ///   This results in a system of two linear equations:<br>
  ///   - (i) @f$ 0 = \vec s(\lambda_0, \mu_0) \cdot \vec e_a = \vec m_ab \cdot
  ///  \vec e_a + \mu_0 \vec e_a \cdot \vec e_b - \lambda_0 @f$ <br>
  ///   - (ii) @f$ 0 = \vec s(\lambda_0, \mu_0) \cdot \vec e_b = \vec m_ab \cdot
  ///  \vec e_b + \mu_0  - \lambda_0 \vec e_b \cdot \vec e_a @f$ <br>
  /// 
  ///   Solving (i), (ii) for @f$ \lambda_0 @f$ and @f$ \mu_0 @f$ yields:
  ///   - @f$ \lambda_0 = \frac{(\vec m_ab \cdot \vec e_a)-(\vec m_ab \cdot \vec
  ///  e_b)(\vec e_a \cdot \vec e_b)}{1-(\vec e_a \cdot \vec e_b)^2} @f$ <br>
  ///   - @f$ \mu_0 = - \frac{(\vec m_ab \cdot \vec e_b)-(\vec m_ab \cdot \vec
  ///  e_a)(\vec e_a \cdot \vec e_b)}{1-(\vec e_a \cdot \vec e_b)^2} @f$ <br>
   virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      dir,
                       bool                 forceDir,
                       const BoundaryCheck& bchk = true) const override;

  /// the pathCorrection for derived classes with thickness
  /// is by definition 1 for StraightLineSurfaces                     
  virtual double
  pathCorrection(const Vector3D&, const Vector3D&) const override
  {
    return 1.;
  }

  /// This method checks if the provided GlobalPosition is inside the assigned
  ///  straw radius, but
  ///  no check is done whether the GlobalPosition is inside bounds or not.
  ///  It overwrites isOnSurface from Base Class as it saves the time of sign
  ///  determination.
  virtual bool
  isOnSurface(const Vector3D&      gpos,
              const BoundaryCheck& bchk = true) const override;

  ///This method returns the bounds of the Surface by reference */
  virtual const SurfaceBounds&
  bounds() const override;

  /// Return properly formatted class name for screen output */
  virtual std::string
  name() const override
  {
    return "Acts::StraightLineSurface";
  };

protected:  
  std::shared_ptr<const CylinderBounds> m_bounds;         ///< bounds (shared)
};

inline StraightLineSurface*
StraightLineSurface::clone(const Transform3D* shift) const
{
  if (shift) new StraightLineSurface(*this, *shift);
  return new StraightLineSurface(*this);
}

inline const SurfaceBounds&
StraightLineSurface::bounds() const
{
  return (*m_bounds.get());
}

inline const Vector3D
StraightLineSurface::lineDirection() const
{
  return std::move(Vector3D(transform().rotation().col(2)));
}

}  // end of namespace

#endif  // ACTS_SURFACESSTRAIGHTLINESURFACE_H
