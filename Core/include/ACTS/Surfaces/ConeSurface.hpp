// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_CONESURFACE_H
#define ACTS_SURFACES_CONESURFACE_H 1

#include "ACTS/Surfaces/ConeBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {

/// @class ConeSurface
///
/// Class for a conical surface in the Tracking geometry.
/// It inherits from Surface.
///
/// The ConeSurface is special since no corresponding
/// Track parameters exist since they're numerical instable
/// at the tip of the cone.
/// Propagations to a cone surface will be returned in
/// curvilinear coordinates.

class ConeSurface : public Surface
{
public:
  /// Default constructor is deleted
  ConeSurface() = delete;

  /// Constructor form HepTransform and an opening angle
  ///
  /// @param htrans is the transform to place to cone in a 3D frame
  /// @param alpha is the opening angle of the cone
  /// @param symmetric indicates if the cones are built to +/1 z
  ConeSurface(std::shared_ptr<Transform3D> htrans,
              double                       alpha,
              bool                         symmetric = false);

  /// Constructor form HepTransform and an opening angle
  ///
  /// @param htrans is the transform that places the cone in the global frame
  /// @param alpha is the opening angle of the cone
  /// @param locZmin is the z range over which the cone spans
  /// @param locZmax is the z range over which the cone spans
  /// @param halfPhi is the openen angle for cone ssectors
  ConeSurface(std::shared_ptr<Transform3D> htrans,
              double                       alpha,
              double                       locZmin,
              double                       locZmax,
              double                       halfPhi = M_PI);

  /// Constructor from HepTransform and ConeBounds
  ///
  /// @param htrans is the transform that places the cone in the global frame
  /// @param cbounds is the boundary class, the bounds must exit
  ConeSurface(std::shared_ptr<Transform3D>      htrans,
              std::shared_ptr<const ConeBounds> cbounds);

  /// Copy constructor
  ///
  /// @param csf is the source cone surface
  ConeSurface(const ConeSurface& csf);

  /// Copy constructor - with shift
  ///
  /// @param csf is the source cone surface
  /// @param htrans is the additional transfrom applied after copying
  ConeSurface(const ConeSurface& csf, const Transform3D& htrans);

  /// Destructor
  virtual ~ConeSurface();

  /// Assignment operator
  ///
  /// @param csf is the source surface for the assignment
  ConeSurface&
  operator=(const ConeSurface& csf);

  /// Implicit Constructor
  ///
  /// @param shift is the optional shift applied after cloning
  virtual ConeSurface*
  clone(const Transform3D* shift = nullptr) const final override;

  /// The binning position method - is overloaded for r-type binning
  ///
  /// @param bValue defines the type of binning to be applied in the global
  /// frame
  ///
  /// @return The return type is a vector for positioning in the global frame
  virtual const Vector3D
  binningPosition(BinningValue bValue) const final override;

  /// Return the surface type
  virtual SurfaceType
  type() const override
  {
    return Surface::Cone;
  }

  /// Return the measurement frame - this is needed for alignment, in particular
  ///  for StraightLine and Perigee Surface
  ///  - the default implementation is the the RotationMatrix3D of the transform
  ///
  /// @param gpos is the global position where the measurement frame is
  /// constructed
  /// @param mom is the momentum used for the measurement frame construction
  ///
  /// @return matrix that indicates the measurement frame
  const RotationMatrix3D
  measurementFrame(const Vector3D& gpos, const Vector3D& mom) const final;

  /// Return method for surface normal information
  ///
  /// @param lpos is the local position on the cone for which the normal vector
  /// is requested
  ///
  /// @return Vector3D normal vector in global frame
  const Vector3D
  normal(const Vector2D& lpos) const final;

  /// Return method for surface normal information
  ///
  /// @param gpos is the global position on the cone for which the normal vector
  /// is requested
  ///
  /// @return Vector3D normal vector in global frame
  const Vector3D
  normal(const Vector3D& gpos) const final;

  // Return method for the rotational symmetry axis
  ///
  // @return This returns the local z axis
  virtual const Vector3D
  rotSymmetryAxis() const;

  /// This method returns the ConeBounds by reference
  virtual const ConeBounds&
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
  ///
  /// @return is a boolean indicating if the transformation succeeded
  virtual bool
  globalToLocal(const Vector3D& gpos,
                const Vector3D& mom,
                Vector2D&       lpos) const final override;

  /// straight line intersection schema - provides closest intersection and
  /// (signed) path length
  ///
  /// @param gpos is the start position for the intersection
  /// @param dir is the start direction for the intersection
  /// @param forceDir is the flag to force to go along the forward direction
  /// @param bcheck is the boundary check to be used in this directive
  ///
  /// <b>mathematical motivation:</b>
  ///
  ///   The calculation will be done in the 3-dim frame of the cone,
  ///   i.e. the symmetry axis of the cone is the z-axis, x- and y-axis are
  /// perpenticular
  ///   to the the z-axis. In this frame the cone is centered around the origin.
  ///   Therefore the two points describing the line have to be first
  ///   recalculated
  /// into the new frame.
  ///   Suppose, this is done, the points of intersection can be
  ///   obtained as follows:<br>
  ///
  ///   The cone is described by the implicit equation
  ///   @f$x^2 + y^2 = z^2 \tan \alpha@f$
  ///   where @f$\alpha@f$ is opening half-angle of the cone  the and
  ///   the line by the parameter equation (with @f$t@f$ the
  ///   parameter and @f$x_1@f$ and @f$x_2@f$ are points on the line)
  ///   @f$(x,y,z) = \vec x_1 + (\vec x_2 - \vec x_2) t @f$.
  ///   The intersection is the given to the value of @f$t@f$ where
  ///   the @f$(x,y,z)@f$ coordinates of the line satisfy the implicit
  ///   equation of the cone. Inserting the expression for the points
  ///   on the line into the equation of the cone and rearranging to
  ///   the form of a  gives (letting @f$ \vec x_d = \frac{\vec x_2 - \vec
  ///   x_1}{|\vec x_2 - \vec x_1|} @f$):
  ///   @f$t^2 (x_d^2 + y_d^2 - z_d^2 \tan^2 \alpha) + 2 t (x_1 x_d +
  ///   y_1 y_d - z_1 z_d \tan^2 \alpha) + (x_1^2 + y_1^2 - z_1^2
  ///   \tan^2 \alpha) = 0 @f$
  ///   Solving the above for @f$t@f$ and putting the values into the
  ///   equation of the line gives the points of intersection. @f$t@f$
  ///   is also the length of the path, since we normalized @f$x_d@f$
  ///   to be unit length.
  ///
  /// @return is the Intersection object
  virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      dir,
                       bool                 forceDir = false,
                       const BoundaryCheck& bcheck   = false) const final override;

  /// the pathCorrection for derived classes with thickness
  ///
  /// @param gpos is the global potion at the correction point
  /// @param mom is the momentum at the correction point
  ///
  /// @return is the path correction due to incident angle
  virtual double
  pathCorrection(const Vector3D& gpos, const Vector3D& mom) const final override;

  /// Return properly formatted class name for screen output
  virtual std::string
  name() const override
  {
    return "Acts::ConeSurface";
  }

protected:
  std::shared_ptr<const ConeBounds> m_bounds;  ///< bounds (shared)
};

inline ConeSurface*
ConeSurface::clone(const Transform3D* shift) const
{
  if (shift) new ConeSurface(*this, *shift);
  return new ConeSurface(*this);
}

inline const Vector3D
ConeSurface::normal(const Vector2D& lp) const
{
  // (cos phi cos alpha, sin phi cos alpha, sgn z sin alpha)
  double phi = lp[Acts::eLOC_RPHI] / (bounds().r(lp[Acts::eLOC_Z])),
         sgn = lp[Acts::eLOC_Z] > 0 ? -1. : +1.;
  Vector3D localNormal(cos(phi) * bounds().cosAlpha(),
                       sin(phi) * bounds().cosAlpha(),
                       sgn * bounds().sinAlpha());
  return m_transform ? Vector3D(transform().linear() * localNormal)
                     : localNormal;
}

inline const Vector3D
ConeSurface::normal(const Vector3D& gpos) const
{
  // get it into the cylinder frame if needed
  // @todo respect opening angle
  Vector3D pos3D = gpos;
  if (m_transform || m_associatedDetElement) {
    pos3D     = transform().inverse() * gpos;
    pos3D.z() = 0;
  }
  return pos3D.unit();
}

inline const ConeBounds&
ConeSurface::bounds() const
{
  // is safe because no constructor w/o bounds exists
  return (*m_bounds.get());
}

}  // end of namespace

#endif  // ACTS_SURFACESCONESURFACE_H
