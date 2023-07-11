// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/PointBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <string>

namespace Acts {

class PointBounds;
class DetectorElementBase;
class SurfaceBounds;

/// @class PointSurface
///
/// Base class for a point surfaces in the TrackingGeometry to allow propagation
/// to the closest point of approach of a 3D point. It inherits from Surface.
///
/// @note It leaves the type() method virtual, so it can not be instantiated
class PointSurface : public Surface {
#ifndef DOXYGEN
  friend Surface;
#endif

 protected:
  /// Constructor from Translation3 and bounds
  ///
  /// @param transform The transform that positions the surface in the global
  /// frame
  /// @param radius The straw radius
  /// @param halez The half length in z
  PointSurface(const Translation3& transform, double radius);

  /// Constructor from Translation3 and a shared bounds object
  ///
  /// @param transform The transform that positions the surface in the global
  /// frame
  /// @param pbounds The bounds describing the straw dimensions, can be
  /// optionally nullptr
  PointSurface(const Translation3& transform,
               std::shared_ptr<const PointBounds> pbounds = nullptr);

  /// Constructor from DetectorElementBase : Element proxy
  ///
  /// @param pbounds The bounds describing the straw dimensions
  /// @param detelement for which this surface is (at least) one representation
  PointSurface(std::shared_ptr<const PointBounds> pbounds,
               const DetectorElementBase& detelement);

  /// Copy constructor
  ///
  /// @param other The source surface for copying
  PointSurface(const PointSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source cone surface
  /// @param shift is the additional transform applied after copying
  PointSurface(const GeometryContext& gctx, const PointSurface& other,
               const Transform3& shift);

 public:
  ~PointSurface() override = default;
  PointSurface() = delete;

  /// Assignment operator
  ///
  /// @param other is the source surface dor copying
  PointSurface& operator=(const PointSurface& other);

  /// Normal vector return
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position is ignored
  ///
  /// @return a Vector3 by value
  Vector3 normal(const GeometryContext& gctx,
                 const Vector2& lposition) const final;

  /// Normal vector return without argument
  using Surface::normal;

  /// The binning position is the position calculated
  /// for a certain binning type
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the binning type to be used
  ///
  /// @return position that can beused for this binning
  Vector3 binningPosition(const GeometryContext& gctx,
                          BinningValue bValue) const final;

  /// Return the measurement frame
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position where the measurement frame is
  /// constructed
  /// @param direction is the momentum direction used for the measurement frame
  /// construction
  ///
  /// @return is a rotation matrix that indicates the measurement frame
  RotationMatrix3 referenceFrame(const GeometryContext& gctx,
                                 const Vector3& position,
                                 const Vector3& direction) const final;

  /// Calculate the jacobian from local to global which the surface knows best,
  /// hence the calculation is done here.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param boundParams is the bound parameters vector
  ///
  /// @return Jacobian from local to global
  BoundToFreeMatrix boundToFreeJacobian(
      const GeometryContext& gctx, const BoundVector& boundParams) const final;

  /// Calculate the derivative of path length at the geometry constraint or
  /// point-of-closest-approach w.r.t. free parameters
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param parameters is the free parameters
  ///
  /// @return Derivative of path length w.r.t. free parameters
  FreeToPathMatrix freeToPathDerivative(
      const GeometryContext& gctx, const FreeVector& parameters) const final;

  /// Local to global transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position to be transformed
  /// @param direction is the global momentum direction (used to sign the closest approach)
  ///
  /// @return global position by value
  Vector3 localToGlobal(const GeometryContext& gctx, const Vector2& lposition,
                        const Vector3& direction) const final;

  /// This method is the true global->local transformation.
  ///
  /// Makes use of globalToLocal and indicates the sign of the Acts::eBoundLoc0
  /// by the given momentum direction.
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param direction global 3D momentum direction (optionally ignored)
  /// @param tolerance (unused)
  ///
  /// @return a Result<Vector2> which can be !ok() if the operation fails
  Result<Vector2> globalToLocal(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction,
      double tolerance = s_onSurfaceTolerance) const final;

  /// @brief Straight line intersection schema
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The global position as a starting point
  /// @param direction The global direction at the starting point
  ///        @note expected to be normalized
  /// @param bcheck The boundary check directive for the estimate
  /// @param tolerance the tolerance used for the intersection
  ///
  /// @return is the intersection object
  SurfaceIntersection intersect(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction, const BoundaryCheck& bcheck = false,
      ActsScalar tolerance = s_onSurfaceTolerance) const final;

  /// The pathCorrection for derived classes with thickness is by definition 1
  /// for PointSurface.
  ///
  /// @note input parameters are ignored
  /// @note there's no material associated to the point surface
  double pathCorrection(const GeometryContext& gctx, const Vector3& position,
                        const Vector3& direction) const override;

  /// This method returns the bounds of the Surface by reference */
  const SurfaceBounds& bounds() const final;

  /// Return properly formatted class name for screen output */
  std::string name() const override;

  /// Calculate the derivative of path length at the geometry constraint or
  /// point-of-closest-approach w.r.t. alignment parameters of the surface (i.e.
  /// local frame origin in global 3D Cartesian coordinates and its rotation
  /// represented with extrinsic Euler angles)
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param parameters is the free parameters
  ///
  /// @return Derivative of path length w.r.t. the alignment parameters
  AlignmentToPathMatrix alignmentToPathDerivative(
      const GeometryContext& gctx, const FreeVector& parameters) const final;

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
  std::shared_ptr<const PointBounds> m_bounds;  ///< bounds (shared)

 private:
  /// helper function to apply the globalToLocal with out transform
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position
  /// @param direction is the momentum direction
  /// @param lposition is the local position to be filled
  bool globalToLocalPlain(const GeometryContext& gctx, const Vector3& position,
                          const Vector3& direction, Vector2& lposition) const;
};

}  // namespace Acts
