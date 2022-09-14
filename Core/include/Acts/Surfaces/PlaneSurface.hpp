// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <limits>

namespace Acts {

class DetectorElementBase;

/// @class PlaneSurface
///
/// Class for a planaer in the TrackingGeometry.
///
/// The PlaneSurface extends the Surface class with the possibility to
/// convert local to global positions (vice versa).
///
/// @image html figures/PlaneSurface.png
///
class PlaneSurface : public Surface {
#ifndef DOXYGEN
  friend Surface;
#endif

 protected:
  /// Copy Constructor
  ///
  /// @param other is the source surface for the copy
  PlaneSurface(const PlaneSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source cone surface
  /// @param transf is the additional transfrom applied after copying
  PlaneSurface(const GeometryContext& gctx, const PlaneSurface& other,
               const Transform3& transf);

  /// Dedicated Constructor with normal vector
  /// This is for curvilinear surfaces which are by definition boundless
  ///
  /// @param center is the center position of the surface
  /// @param normal is thenormal vector of the plane surface
  PlaneSurface(const Vector3& center, const Vector3& normal);

  /// Constructor from DetectorElementBase : Element proxy
  ///
  /// @param pbounds are the provided planar bounds (shared)
  /// @param detelement is the linked detector element to this surface
  PlaneSurface(const std::shared_ptr<const PlanarBounds>& pbounds,
               const DetectorElementBase& detelement);

  /// Constructor for Planes with (optional) shared bounds object
  ///
  /// @param htrans transform in 3D that positions this surface
  /// @param pbounds bounds object to describe the actual surface area
  PlaneSurface(const Transform3& htrans,
               std::shared_ptr<const PlanarBounds> pbounds = nullptr);

 public:
  ~PlaneSurface() override = default;
  PlaneSurface() = delete;

  /// Assignment operator
  ///
  /// @param other The source PlaneSurface for assignment
  PlaneSurface& operator=(const PlaneSurface& other);

  /// Normal vector return
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition is the local position is ignored
  ///
  /// return a Vector3 by value
  Vector3 normal(const GeometryContext& gctx,
                 const Vector2& lposition) const final;

  /// Normal vector return without argument
  using Surface::normal;

  /// The binning position is the position calcualted
  /// for a certain binning type
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue is the binning type to be used
  ///
  /// @return position that can beused for this binning
  Vector3 binningPosition(const GeometryContext& gctx,
                          BinningValue bValue) const final;

  /// Return the surface type
  SurfaceType type() const override;

  /// Return method for bounds object of this surfrace
  const SurfaceBounds& bounds() const override;

  /// Local to global transformation
  /// For planar surfaces the momentum is ignroed in the local to global
  /// transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lposition local 2D position in specialized surface frame
  /// @param momentum global 3D momentum representation (optionally ignored)
  ///
  /// @return the global position by value
  Vector3 localToGlobal(const GeometryContext& gctx, const Vector2& lposition,
                        const Vector3& momentum) const override;

  /// Global to local transformation
  /// For planar surfaces the momentum is ignroed in the global to local
  /// transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param momentum global 3D momentum representation (optionally ignored)
  /// method symmetry)
  /// @param tolerance optional tolerance within which a point is considered
  /// valid on surface
  ///
  /// @return a Result<Vector2> which can be !ok() if the operation fails
  Result<Vector2> globalToLocal(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& momentum,
      double tolerance = s_onSurfaceTolerance) const override;

  /// Method that calculates the correction due to incident angle
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param direction global 3D momentum direction (ignored for PlaneSurface)
  /// @note this is the final implementation of the pathCorrection function
  ///
  /// @return a double representing the scaling factor
  double pathCorrection(const GeometryContext& gctx, const Vector3& position,
                        const Vector3& direction) const final;

  /// @brief Straight line intersection
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The start position of the intersection attempt
  /// @param direction The direction of the interesection attempt,
  /// (@note expected to be normalized)
  /// @param bcheck The boundary check directive
  ///
  /// <b>mathematical motivation:</b>
  ///
  /// the equation of the plane is given by: <br>
  /// @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
  /// where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of
  /// the plane,  @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point
  /// on the plane and @f$ \vec x = (x,y,z) @f$ all possible points
  /// on the plane.<br>
  ///
  /// Given a line with:<br>
  /// @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
  /// the solution for @f$ u @f$ can be written:
  /// @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
  /// If the denominator is 0 then the line lies:
  /// - either in the plane
  /// - perpendicular to the normal of the plane
  ///
  /// @return the SurfaceIntersection object
  SurfaceIntersection intersect(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction,
      const BoundaryCheck& bcheck = false) const final;

  /// Return a Polyhedron for the surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lseg Number of segments along curved lines, it represents
  /// the full 2*M_PI coverange, if lseg is set to 1 only the extrema
  /// are given
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& gctx,
                                      size_t lseg) const override;

  /// Return properly formatted class name for screen output
  std::string name() const override;

  /// Calculate the derivative of bound track parameters local position w.r.t.
  /// position in local 3D Cartesian coordinates
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position of the paramters in global
  ///
  /// @return Derivative of bound local position w.r.t. position in local 3D
  /// cartesian coordinates
  ActsMatrix<2, 3> localCartesianToBoundLocalDerivative(
      const GeometryContext& gctx, const Vector3& position) const final;

 protected:
  /// the bounds of this surface
  std::shared_ptr<const PlanarBounds> m_bounds;

 private:
};

}  // end of namespace Acts
