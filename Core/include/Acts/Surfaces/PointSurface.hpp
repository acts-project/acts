// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/PointBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceConcept.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Result.hpp"

#include <iosfwd>
#include <memory>
#include <string>

namespace Acts {

class PointBounds;
class SurfaceBounds;

/// @class PointSurface
///
/// Surface describing the point of closest approach (PCA) of a track to a
/// single fixed 3D point (the surface center). It is the point-analogue of the
/// LineSurface (which represents the PCA to a line / Perigee).
///
/// Geometrically, a PointSurface behaves like a curvilinear plane anchored at
/// the fixed center: a plane through the center whose normal always equals the
/// track momentum direction. The track's crossing of that plane is therefore
/// the PCA to the point. Consequently there is no parallel-line degeneracy: a
/// unique intersection always exists.
///
/// The local coordinates are Cartesian (x, y) in the measurement plane
/// perpendicular to the track direction. Cartesian (not polar r-phi) is used on
/// purpose, because at the point itself (r = 0) the azimuth phi is degenerate.
///
/// The measurement frame is oriented from the track direction like the
/// CurvilinearSurface: the third axis is the track direction, the first axis is
/// built from global Z (with a fall-back to global X when the direction is
/// (anti-)parallel to Z).
///
/// @note A point has no intrinsic orientation, so only translation (center)
/// alignment is supported; rotational alignment derivatives are zero.
class PointSurface : public Surface {
  friend class Surface;

 protected:
  /// Constructor from a global position (unbounded point surface)
  ///
  /// @param center position of the point in the global frame
  explicit PointSurface(const Vector3& center);

  /// Constructor from a global position and a maximum distance
  ///
  /// @param center position of the point in the global frame
  /// @param maxDistance maximum distance from the point (radius of the
  /// PointBounds disc in the measurement plane)
  PointSurface(const Vector3& center, double maxDistance);

  /// Constructor from a Transform3 and optional PointBounds
  ///
  /// @param transform the transform whose translation positions the point
  /// @param pbounds the bounds describing the maximum distance, may be nullptr
  /// for an unbounded point surface
  explicit PointSurface(const Transform3& transform,
                        std::shared_ptr<const PointBounds> pbounds = nullptr);

  /// Copy constructor
  ///
  /// @param other The source surface for copying
  PointSurface(const PointSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source surface
  /// @param shift is the additional transform applied after copying
  PointSurface(const GeometryContext& gctx, const PointSurface& other,
               const Transform3& shift);

 public:
  /// Assignment operator
  ///
  /// @param other is the source surface for copying
  /// @return Reference to this PointSurface after assignment
  PointSurface& operator=(const PointSurface& other);

  /// Return the surface type
  /// @return Surface type identifier for point surfaces
  SurfaceType type() const final { return Surface::Point; }

  /// Return the surface normal, which for a point surface is the momentum
  /// direction (the measurement plane is perpendicular to the direction).
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param pos is the global position (ignored)
  /// @param direction is the global momentum direction
  /// @return the normal vector (equal to @p direction)
  Vector3 normal(const GeometryContext& gctx, const Vector3& pos,
                 const Vector3& direction) const final;

  /// The binning position is the center of the point surface.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param aDir is the axis direction for the reference position request
  /// @return the center position
  Vector3 referencePosition(const GeometryContext& gctx,
                            AxisDirection aDir) const final;

  /// Return the measurement frame, built from the momentum direction like a
  /// curvilinear surface. Columns are (U, V, T) with T the direction.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position (ignored)
  /// @param direction is the momentum direction used to build the frame
  /// @return a rotation matrix that indicates the measurement frame
  RotationMatrix3 referenceFrame(const GeometryContext& gctx,
                                 const Vector3& position,
                                 const Vector3& direction) const final;

  /// Calculate the jacobian from local to global.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position
  /// @param direction global 3D momentum direction
  /// @return Jacobian from local to global
  BoundToFreeMatrix boundToFreeJacobian(const GeometryContext& gctx,
                                        const Vector3& position,
                                        const Vector3& direction) const final;

  /// Calculate the derivative of path length at the point-of-closest-approach
  /// w.r.t. free parameters.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position
  /// @param direction global 3D momentum direction
  /// @return Derivative of path length w.r.t. free parameters
  FreeToPathMatrix freeToPathDerivative(const GeometryContext& gctx,
                                        const Vector3& position,
                                        const Vector3& direction) const final;

  /// Local to global transformation
  ///
  /// @note for point surfaces the momentum direction is used in order to
  /// build the measurement plane
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local (x, y) position to be transformed
  /// @param direction is the global momentum direction
  /// @return global position by value
  Vector3 localToGlobal(const GeometryContext& gctx, const Vector2& lposition,
                        const Vector3& direction) const final;

  /// Global to local transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param direction global 3D momentum direction
  /// @param tolerance the tolerance for the on-surface check
  /// @return A Result<Vector2> which is !ok() if @p position is not the point
  /// of closest approach to the point surface.
  Result<Vector2> globalToLocal(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction,
      double tolerance = s_onSurfaceTolerance) const final;

  /// Calculate the straight-line intersection with the point surface, i.e. the
  /// point of closest approach of the track to the point.
  ///
  /// Given the track (@p position @f$ \vec m @f$, @p direction @f$ \vec e @f$,
  /// normalized), the PCA path length is @f$ u = (\vec c - \vec m) \cdot \vec e
  /// @f$ where @f$ \vec c @f$ is the point, and the intersection point is @f$
  /// \vec m + u \vec e @f$. Unlike the LineSurface there is no parallel
  /// degeneracy.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The global position as a starting point
  /// @param direction The global direction at the starting point
  ///        @note expected to be normalized
  /// @param boundaryTolerance The boundary check directive for the estimate
  /// @param tolerance the tolerance used for the intersection
  /// @return is the intersection object
  MultiIntersection3D intersect(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction,
      const BoundaryTolerance& boundaryTolerance =
          BoundaryTolerance::Infinite(),
      double tolerance = s_onSurfaceTolerance) const final;

  /// The pathCorrection is by definition 1 for point surfaces
  ///
  /// @param gctx Geometry context (ignored)
  /// @param position Position parameter (ignored)
  /// @param direction Direction parameter (ignored)
  /// @return Always returns 1.0 for point surfaces
  double pathCorrection(const GeometryContext& gctx, const Vector3& position,
                        const Vector3& direction) const final;

  /// This method returns the bounds of the surface by reference
  /// @return Reference to the surface bounds
  const SurfaceBounds& bounds() const final;
  /// This method returns the shared_ptr to the PointBounds
  /// @return Shared pointer to the point bounds (may be nullptr)
  const std::shared_ptr<const PointBounds>& boundsPtr() const;
  /// Overwrite the existing surface bounds with new ones
  /// @param newBounds Pointer to the new bounds
  void assignSurfaceBounds(std::shared_ptr<const PointBounds> newBounds);

  /// Return properly formatted class name for screen output
  /// @return String representation of the class name
  std::string name() const final;

  /// Return a Polyhedron for the surface (a non-physical marker at the center)
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param quarterSegments is an ignored parameter
  /// @return A list of vertices and a face/facet description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& gctx,
                                      unsigned int quarterSegments) const final;

  /// Calculate the derivative of path length at the point-of-closest-approach
  /// w.r.t. alignment parameters of the surface. Only the center (translation)
  /// contributes; a point has no orientation so the rotation block is zero.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position
  /// @param direction global 3D momentum direction
  /// @return Derivative of path length w.r.t. the alignment parameters
  AlignmentToPathMatrix alignmentToPathDerivative(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction) const final;

  /// Calculate the derivative of bound track parameters local position w.r.t.
  /// position in local 3D Cartesian coordinates
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position of the parameters in global
  /// @return Derivative of bound local position w.r.t. position in local 3D
  /// cartesian coordinates
  Matrix<2, 3> localCartesianToBoundLocalDerivative(
      const GeometryContext& gctx, const Vector3& position) const final;

 protected:
  std::shared_ptr<const PointBounds>
      m_bounds;  ///< bounds (shared, may be null)

  /// @copydoc Surface::localAxes
  std::array<AxisDirection, 2> localAxes() const final {
    return {AxisDirection::AxisX, AxisDirection::AxisY};
  }

  /// Output Method for std::ostream
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param sl is the ostream to be dumped into
  /// @return ostream object which was streamed into
  std::ostream& toStreamImpl(const GeometryContext& gctx,
                             std::ostream& sl) const final;
};

static_assert(SurfaceConcept<PointSurface>,
              "PointSurface does not fulfill SurfaceConcept");

}  // namespace Acts
