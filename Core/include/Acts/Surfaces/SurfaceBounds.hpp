// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/OstreamFormatter.hpp"

#include <cmath>
#include <ostream>

namespace Acts {

/// @class SurfaceBounds
///
/// Interface for surface bounds.
///
/// Surface bounds provide:
/// - inside() checks
/// - distance to boundary calculations
/// - the BoundsType and a set of parameters to simplify persistency
///
class SurfaceBounds {
 public:
  /// @enum BoundsType
  /// This is nested to the SurfaceBounds, as also VolumeBounds will have
  /// Bounds Type.
  enum class BoundsType : int {
    eCone [[deprecated("Use the enum ::Cone instead")]] = 0,
    eCylinder [[deprecated("Use the enum ::Cylinder instead")]] = 1,
    eDiamond [[deprecated("Use the enum ::Diamond instead")]] = 2,
    eDisc [[deprecated("Use the enum ::Disc instead")]] = 3,
    eEllipse [[deprecated("Use the enum ::Ellipse instead")]] = 4,
    eLine [[deprecated("Use the enum ::Line instead")]] = 5,
    eRectangle [[deprecated("Use the enum ::Rectangle instead")]] = 6,
    eTrapezoid [[deprecated("Use the enum ::Trapezoid instead")]] = 7,
    eTriangle [[deprecated("Use the enum ::Triangle instead")]] = 8,
    eDiscTrapezoid = 9,
    eConvexPolygon [[deprecated("Use the enum ::ConvexPolygon instead")]] = 10,
    eAnnulus [[deprecated("Use the enum ::Annulus instead")]] = 11,
    eBoundless [[deprecated("Use the enum ::Boundless instead")]] = 12,
    eOther [[deprecated("Use the enum ::Other instead")]] = 13,

    Cone = 0,
    Cylinder = 1,
    Diamond = 2,
    Disc = 3,
    Ellipse = 4,
    Line = 5,
    Rectangle = 6,
    Trapezoid = 7,
    Triangle = 8,
    DiscTrapezoid = 9,
    ConvexPolygon = 10,
    Annulus = 11,
    Boundless = 12,
    Other = 13
  };

  using enum BoundsType;

  virtual ~SurfaceBounds() = default;

  /// Return the bounds type - for persistency optimization
  /// @return the bounds type
  virtual BoundsType type() const = 0;

  /// Check if the bound coordinates are cartesian
  /// @return true if the bound coordinates are cartesian
  virtual bool isCartesian() const = 0;

  /// Computes the bound to cartesian jacobian at a given local position
  /// @param lposition is the local position at which the jacobian is computed
  /// @return the bound to cartesian jacobian
  virtual SquareMatrix2 boundToCartesianJacobian(
      const Vector2& lposition) const = 0;

  /// Computes the bound to cartesian metric at a given local position
  /// @param lposition is the local position at which the metric is computed
  /// @return the bound to cartesian metric
  virtual SquareMatrix2 boundToCartesianMetric(
      const Vector2& lposition) const = 0;

  /// Access method for bound values, this is a dynamically sized
  /// vector containing the parameters needed to describe these bounds
  /// @return of the stored values for this SurfaceBounds object
  virtual std::vector<double> values() const = 0;

  /// Inside check for the bounds object
  /// @param lposition is the local position
  /// @return true if the local position is inside the bounds
  virtual bool inside(const Vector2& lposition) const = 0;

  /// Calculates the closest point on the bounds to a given local position
  /// @param lposition is the local position
  /// @param metric to be used for the distance calculation
  /// @return the closest point on the bounds
  virtual Vector2 closestPoint(const Vector2& lposition,
                               const SquareMatrix2& metric) const = 0;

  /// Calculates the distance to the bounds from a given local position
  /// @param lposition is the local position
  /// @return the distance to the bounds
  virtual double distance(const Vector2& lposition) const {
    SquareMatrix2 metric = boundToCartesianMetric(lposition);

    Vector2 closest = closestPoint(lposition, metric);
    Vector2 diff = closest - lposition;
    return std::sqrt((diff.transpose() * metric * diff)(0, 0));
  }

  /// Inside check for the bounds object given a boundary tolerance.
  /// @param lposition is the local position
  /// @param boundaryTolerance is the boundary tolerance object
  /// @return true if the local position is inside the bounds and tolerance
  virtual bool inside(const Vector2& lposition,
                      const BoundaryTolerance& boundaryTolerance) const;

  /// Calculate the center of the surface bounds in local coordinates
  ///
  /// This method returns a representative center point of the bounds region.
  /// The exact definition varies by bounds type and coordinate system:
  ///
  /// **Cartesian bounds** (Rectangle, Diamond, Trapezoid):
  /// - Returns the geometric center or center of symmetry
  /// - For symmetric shapes: center of bounding box or origin (0,0)
  ///
  /// **Polar/Cylindrical bounds** (Radial, Cylinder, Cone):
  /// - Returns (r, phi) where r is average radius, phi is average angle
  /// - Coordinates are in the bounds' natural coordinate system
  ///
  /// **Complex bounds** (Annulus, ConvexPolygon):
  /// - Annulus: Pre-calculated from corner vertices (accounts for coordinate
  /// transforms)
  /// - Polygon: Average of all vertices (vertex centroid, not area centroid)
  ///
  /// **Infinite bounds**: Returns conceptual center at (0,0)
  ///
  /// @note The returned point is guaranteed to be a reasonable representative
  /// center, but may not be the true geometric centroid for all shapes.
  ///
  /// @return Vector2 representing the center position in local coordinates
  virtual Vector2 center() const = 0;

  /// Output Method for std::ostream, to be overloaded by child classes
  /// @param os is the outstream in which the string dump is done
  /// @return Modified ostream for chaining
  virtual std::ostream& toStream(std::ostream& os) const = 0;

  friend bool operator==(const SurfaceBounds& lhs, const SurfaceBounds& rhs) {
    if (&lhs == &rhs) {
      return true;
    }
    return (lhs.type() == rhs.type()) && (lhs.values() == rhs.values());
  }

  friend std::ostream& operator<<(std::ostream& os, const SurfaceBounds& sb) {
    return sb.toStream(os);
  }
};

/// Stream operator for SurfaceBounds::BoundsType
std::ostream& operator<<(std::ostream& os, const SurfaceBounds::BoundsType& bt);

}  // namespace Acts

namespace std {
template <>
struct formatter<Acts::SurfaceBounds::BoundsType>
    : Acts::detail::OstreamFormatter {};
}  // namespace std
