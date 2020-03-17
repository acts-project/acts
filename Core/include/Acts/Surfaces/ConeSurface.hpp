// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"

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

class ConeSurface : public Surface {
  friend Surface;

 protected:
  /// Constructor form HepTransform and an opening angle
  ///
  /// @param htrans is the transform to place to cone in a 3D frame
  /// @param alpha is the opening angle of the cone
  /// @param symmetric indicates if the cones are built to +/1 z
  ConeSurface(std::shared_ptr<const Transform3D> htrans, double alpha,
              bool symmetric = false);

  /// Constructor form HepTransform and an opening angle
  ///
  /// @param htrans is the transform that places the cone in the global frame
  /// @param alpha is the opening angle of the cone
  /// @param zmin is the z range over which the cone spans
  /// @param zmax is the z range over which the cone spans
  /// @param halfPhi is the openen angle for cone ssectors
  ConeSurface(std::shared_ptr<const Transform3D> htrans, double alpha,
              double zmin, double zmax, double halfPhi = M_PI);

  /// Constructor from HepTransform and ConeBounds
  ///
  /// @param htrans is the transform that places the cone in the global frame
  /// @param cbounds is the boundary class, the bounds must exit
  ConeSurface(std::shared_ptr<const Transform3D> htrans,
              const std::shared_ptr<const ConeBounds>& cbounds);

  /// Copy constructor
  ///
  /// @param other is the source cone surface
  ConeSurface(const ConeSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source cone surface
  /// @param transf is the additional transfrom applied after copying
  ConeSurface(const GeometryContext& gctx, const ConeSurface& other,
              const Transform3D& transf);

 public:
  /// Destructor - defaulted
  ~ConeSurface() override = default;

  /// Default Constructor - deleted
  ConeSurface() = delete;

  /// Assignment operator
  ///
  /// @param other is the source surface for the assignment
  ConeSurface& operator=(const ConeSurface& other);

  /// Clone method into a concrete type of ConeSurface with shift
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param shift applied to the surface
  std::shared_ptr<ConeSurface> clone(const GeometryContext& gctx,
                                     const Transform3D& shift) const;

  /// The binning position method - is overloaded for r-type binning
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue defines the type of binning applied in the global frame
  ///
  /// @return The return type is a vector for positioning in the global frame
  const Vector3D binningPosition(const GeometryContext& gctx,
                                 BinningValue bValue) const final;

  /// Return the surface type
  SurfaceType type() const override;

  /// Return the measurement frame - this is needed for alignment, in particular
  ///  for StraightLine and Perigee Surface
  ///  - the default implementation is the the RotationMatrix3D of the transform
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position where the measurement frame is
  /// constructed
  /// @param momentum is the momentum used for the measurement frame
  /// construction
  /// @return matrix that indicates the measurement frame
  const RotationMatrix3D referenceFrame(const GeometryContext& gctx,
                                        const Vector3D& position,
                                        const Vector3D& momentum) const final;

  /// Return method for surface normal information
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position at normal vector request
  /// @return Vector3D normal vector in global frame
  const Vector3D normal(const GeometryContext& gctx,
                        const Vector2D& lposition) const final;

  /// Return method for surface normal information
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position as normal vector base
  /// @return Vector3D normal vector in global frame
  const Vector3D normal(const GeometryContext& gctx,
                        const Vector3D& position) const final;

  /// Normal vector return without argument
  using Surface::normal;

  // Return method for the rotational symmetry axis
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  // @return This returns the local z axis
  virtual const Vector3D rotSymmetryAxis(const GeometryContext& gctx) const;

  /// This method returns the ConeBounds by reference
  const ConeBounds& bounds() const final;

  /// Local to global transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position to be transformed
  /// @param momentum is the global momentum (ignored in this operation)
  /// @param position is the global position which is filled
  void localToGlobal(const GeometryContext& gctx, const Vector2D& lposition,
                     const Vector3D& momentum, Vector3D& position) const final;

  /// Global to local transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position to be transformed
  /// @param momentum is the global momentum (ignored in this operation)
  /// @param lposition is the local position to be filled
  /// @return is a boolean indicating if the transformation succeeded
  bool globalToLocal(const GeometryContext& gctx, const Vector3D& position,
                     const Vector3D& momentum, Vector2D& lposition) const final;

  /// @brief Straight line intersection schema - provides closest intersection
  /// and (signed) path length
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The start position for the intersection
  /// @param direciton The start direction for the intersection (expected
  /// normalized)
  /// @param bcheck The boundary check to be used in this directive
  ///
  /// @return is the Intersection object
  Intersection intersectionEstimate(
      const GeometryContext& gctx, const Vector3D& position,
      const Vector3D& direction,
      const BoundaryCheck& bcheck = false) const final;

  /// Straight line intersection schema from position/direction
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position to start from
  /// @param direction The direction at start
  /// @param bcheck the Boundary Check
  ///
  /// If possible returns both solutions for the cylinder
  ///
  /// @return SurfaceIntersection object (contains intersection & surface)
  SurfaceIntersection intersect(const GeometryContext& gctx,
                                const Vector3D& position,
                                const Vector3D& direction,
                                const BoundaryCheck& bcheck) const final;

  /// The pathCorrection for derived classes with thickness
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global potion at the correction point
  /// @param direction is the momentum direction at the correction point
  /// @return is the path correction due to incident angle
  double pathCorrection(const GeometryContext& gctx, const Vector3D& position,
                        const Vector3D& direction) const final;

  /// Return a Polyhedron for the surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lseg Number of segments along curved lines, it represents
  /// the full 2*M_PI coverange, if lseg is set to 1 only the extrema
  /// are given
  /// @note that a surface transform can invalidate the extrema
  /// in the transformed space
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& gctx,
                                      size_t lseg) const override;

  /// Return properly formatted class name for screen output
  std::string name() const override;

 protected:
  std::shared_ptr<const ConeBounds> m_bounds;  ///< bounds (shared)

 private:
  /// Implementation of the intersection solver
  ///
  /// <b>mathematical motivation:</b>
  ///
  ///   The calculation will be done in the 3-dim frame of the cone,
  ///   i.e. the symmetry axis of the cone is the z-axis, x- and y-axis are
  /// perpendicular
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
  /// @return the quadratic equation
  detail::RealQuadraticEquation intersectionSolver(
      const GeometryContext& gctx, const Vector3D& position,
      const Vector3D& direction) const;

  /// Clone method implementation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param shift applied to the surface
  ConeSurface* clone_impl(const GeometryContext& gctx,
                          const Transform3D& shift) const override;
};

#include "Acts/Surfaces/detail/ConeSurface.ipp"

}  // namespace Acts
