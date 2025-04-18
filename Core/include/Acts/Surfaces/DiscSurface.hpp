// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceConcept.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Concepts.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <string>

namespace Acts {

class DetectorElementBase;
class DiscBounds;
class SurfaceBounds;

/// @class DiscSurface
///
/// Class for a disc surface (or a segment thereof)
///
/// The DiscSurface is defined by the local polar coordinates @f$ (r,phi) @f$.
///
/// The surface transform positions the disc such that the origin
/// is at @f$ r=0 @f$, independent of the provided \c DiscBounds.
/// The normal vector of the disc (i.e., the local @f$z@f$-axis) is given by
/// @f$ \vec e_{z} = \vec e_{r} \times\vec e_{phi} @f$.
///
/// The disc surface The only surface type for which the
/// covariance matrix is NOT given in the reference frame.
/// A conversion from polar to cartesian coordinates needs
/// to happen to transfer the local coordinates onto the
/// cartesian reference frame coordinates.
///
/// @image html DiscSurface.png
///
class DiscSurface : public RegularSurface {
  friend class Surface;

 protected:
  /// Constructor for Discs from Transform3, \f$ r_{min}, r_{max} \f$
  ///
  /// @param transform is transform that places the disc in the global 3D space
  /// @param rmin The inner radius of the disc surface
  /// @param rmax The outer radius of the disc surface
  /// @param hphisec The opening angle of the disc surface and is optional
  ///        the default is a full disc
  DiscSurface(const Transform3& transform, double rmin, double rmax,
              double hphisec = M_PI);

  /// Constructor for Discs from Transform3, \f$ r_{min}, r_{max}, hx_{min},
  /// hx_{max} \f$
  /// This is n this case you have DiscTrapezoidBounds
  ///
  /// @param transform is transform that places the disc in the global 3D space
  /// @param minhalfx The half length in x at minimal r
  /// @param maxhalfx The half length in x at maximal r
  /// @param minR The outer radius of the disc surface
  /// @param maxR The inner radius of the disc surface
  /// @param avephi The position in phi (default is 0.)
  /// @param stereo The optional stereo angle
  DiscSurface(const Transform3& transform, double minhalfx, double maxhalfx,
              double minR, double maxR, double avephi = 0., double stereo = 0.);

  /// Constructor for Discs from Transform3 and shared DiscBounds
  ///
  /// @param transform The transform that positions the disc in global 3D
  /// @param dbounds The disc bounds describing the surface coverage
  DiscSurface(const Transform3& transform,
              std::shared_ptr<const DiscBounds> dbounds = nullptr);

  /// Constructor from DetectorElementBase : Element proxy
  ///
  /// @param dbounds The disc bounds describing the surface coverage
  /// @param detelement The detector element represented by this surface
  DiscSurface(std::shared_ptr<const DiscBounds> dbounds,
              const DetectorElementBase& detelement);

  /// Copy Constructor
  ///
  /// @param other The source surface for the copy
  DiscSurface(const DiscSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source cone surface
  /// @param shift is the additional transform applied after copying
  DiscSurface(const GeometryContext& gctx, const DiscSurface& other,
              const Transform3& shift);

 public:
  ~DiscSurface() override = default;
  DiscSurface() = delete;

  /// Assignment operator
  ///
  /// @param other The source sourface for the assignment
  DiscSurface& operator=(const DiscSurface& other);

  /// Return the surface type
  SurfaceType type() const override;

  // User overloads from `RegularSurface`
  using RegularSurface::globalToLocal;
  using RegularSurface::localToGlobal;
  using RegularSurface::normal;

  /// Normal vector return
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition The local position is ignored
  ///
  /// @return a Vector3 by value
  Vector3 normal(const GeometryContext& gctx,
                 const Vector2& lposition) const final;

  /// Get the normal vector of this surface at a given global position
  /// @note The @p position is required to be on-surface.
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global positiono (for @ref DiscSurface this is ignored)
  /// @return The normal vector
  Vector3 normal(const GeometryContext& gctx,
                 const Vector3& position) const final;

  /// Get the normal vector, independent of the location
  /// @param gctx The current geometry context object, e.g. alignment
  /// @return The normal vector
  Vector3 normal(const GeometryContext& gctx) const;

  /// The binning position The position calculated
  /// for a certain binning type
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue The binning type to be used
  ///
  /// @return position that can beused for this binning
  Vector3 binningPosition(const GeometryContext& gctx,
                          BinningValue bValue) const final;

  /// This method returns the bounds by reference
  const SurfaceBounds& bounds() const final;

  /// Local to global transformation
  /// For planar surfaces the momentum direction is ignored in the local to
  /// global transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition local 2D position in specialized surface frame
  ///
  /// @return global position by value
  Vector3 localToGlobal(const GeometryContext& gctx,
                        const Vector2& lposition) const final;

  /// Global to local transformation
  /// @note the direction is ignored for Disc surfaces in this calculateion
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param tolerance optional tolerance within which a point is considered
  /// valid on surface
  ///
  /// @return a Result<Vector2> which can be !ok() if the operation fails
  Result<Vector2> globalToLocal(
      const GeometryContext& gctx, const Vector3& position,
      double tolerance = s_onSurfaceTolerance) const final;

  /// Special method for DiscSurface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lpolar is a local position in polar coordinates
  ///
  /// @return values is local 2D position in cartesian coordinates  @todo check
  Vector2 localPolarToCartesian(const Vector2& lpolar) const;

  /// Special method for Disc surface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lcart is local 2D position in cartesian coordinates
  ///
  /// @return value is a local position in polar coordinates
  Vector2 localCartesianToPolar(const Vector2& lcart) const;

  /// Special method for DiscSurface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param locpol is a local position in polar coordinates
  ///
  /// @return values is local 2D position in cartesian coordinates
  Vector2 localPolarToLocalCartesian(const Vector2& locpol) const;

  /// Special method for DiscSurface :  local<->global transformation when
  /// provided cartesian coordinates
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is local 2D position in cartesian coordinates
  ///
  /// @return value is a global cartesian 3D position
  Vector3 localCartesianToGlobal(const GeometryContext& gctx,
                                 const Vector2& lposition) const;

  /// Special method for DiscSurface : global<->local from cartesian coordinates
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is a global cartesian 3D position
  /// @param tol The absolute tolerance parameter
  ///
  /// @return value is a local polar
  Vector2 globalToLocalCartesian(const GeometryContext& gctx,
                                 const Vector3& position,
                                 double tol = 0.) const;

  /// Calculate the jacobian from local to global which the surface knows best,
  /// hence the calculation is done here.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position
  /// @param direction global 3D momentum direction
  ///
  /// @return Jacobian from local to global
  BoundToFreeMatrix boundToFreeJacobian(const GeometryContext& gctx,
                                        const Vector3& position,
                                        const Vector3& direction) const final;

  /// Calculate the jacobian from global to local which the surface knows best,
  /// hence the calculation is done here.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position
  /// @param direction global 3D momentum direction
  ///
  /// @return Jacobian from global to local
  FreeToBoundMatrix freeToBoundJacobian(const GeometryContext& gctx,
                                        const Vector3& position,
                                        const Vector3& direction) const final;

  /// Path correction due to incident of the track
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The global position as a starting point
  /// @param direction The global momentum direction at the starting point
  /// @return The correction factor due to incident
  double pathCorrection(const GeometryContext& gctx, const Vector3& position,
                        const Vector3& direction) const final;

  /// @brief Straight line intersection schema
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The global position as a starting point
  /// @param direction The global direction at the starting point
  ///        @note expected to be normalized (no checking)
  /// @param bcheck The boundary check prescription
  /// @param tolerance the tolerance used for the intersection
  ///
  /// <b>Mathematical motivation:</b>
  ///
  /// the equation of the plane is given by: <br>
  /// @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
  /// where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of
  /// the plane, @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point on
  /// the plane and @f$ \vec x = (x,y,z) @f$ all possible points
  /// on the plane.<br>
  /// Given a line with:<br>
  /// @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
  /// the solution for @f$ u @f$ can be written:
  /// @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
  /// If the denominator is 0 then the line lies:
  /// - either in the plane
  /// - perpendicular to the normal of the plane
  ///
  /// @return The @c SurfaceMultiIntersection object
  SurfaceMultiIntersection intersect(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction,
      const BoundaryCheck& bcheck = BoundaryCheck(false),
      ActsScalar tolerance = s_onSurfaceTolerance) const final;

  /// Implement the binningValue
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the dobule in which you want to bin
  ///
  /// @note This calls the parent method except for binR
  ///
  /// @return float to be used for the binning schema
  double binningPositionValue(const GeometryContext& gctx,
                              BinningValue bValue) const final;

  /// Return properly formatted class name for screen output
  std::string name() const override;

  /// Return a Polyhedron for the surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lseg Number of segments along curved lines, it represents
  /// the full 2*M_PI coverange, if lseg is set to 1 only the extrema
  /// are given
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& gctx,
                                      std::size_t lseg) const override;

  /// Calculate the derivative of bound track parameters local position w.r.t.
  /// position in local 3D Cartesian coordinates
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position of the parameters in global
  ///
  /// @return Derivative of bound local position w.r.t. position in local 3D
  /// cartesian coordinates
  ActsMatrix<2, 3> localCartesianToBoundLocalDerivative(
      const GeometryContext& gctx, const Vector3& position) const final;

 protected:
  std::shared_ptr<const DiscBounds> m_bounds;  ///< bounds (shared)
};

ACTS_STATIC_CHECK_CONCEPT(RegularSurfaceConcept, DiscSurface);

}  // namespace Acts
