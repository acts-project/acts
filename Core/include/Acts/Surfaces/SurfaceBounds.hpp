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

#include <cmath>
#include <optional>
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
  enum BoundsType : int {
    eCone = 0,
    eCylinder = 1,
    eDiamond = 2,
    eDisc = 3,
    eEllipse = 4,
    eLine = 5,
    eRectangle = 6,
    eTrapezoid = 7,
    eTriangle = 8,
    eDiscTrapezoid = 9,
    eConvexPolygon = 10,
    eAnnulus = 11,
    eBoundless = 12,
    eOther = 13
  };

  virtual ~SurfaceBounds() = default;

  /// Return the bounds type - for persistency optimization
  ///
  /// @return is a BoundsType enum
  virtual BoundsType type() const = 0;

  virtual bool isCartesian() const = 0;

  virtual SquareMatrix2 boundToCartesianJacobian(
      const Vector2& lposition) const = 0;

  virtual SquareMatrix2 cartesianToBoundJacobian(
      const Vector2& lposition) const = 0;

  virtual SquareMatrix2 boundToCartesianMetric(
      const Vector2& lposition) const = 0;

  /// Access method for bound values, this is a dynamically sized
  /// vector containing the parameters needed to describe these bounds
  ///
  /// @return of the stored values for this SurfaceBounds object
  virtual std::vector<double> values() const = 0;

  virtual bool inside(const Vector2& lposition) const = 0;

  virtual Vector2 closestPoint(
      const Vector2& lposition,
      const std::optional<SquareMatrix2>& metric) const = 0;

  virtual double distance(const Vector2& lposition,
                          const std::optional<SquareMatrix2>& metric) const {
    Vector2 closest = closestPoint(lposition, metric);
    Vector2 diff = closest - lposition;
    return std::sqrt((diff.transpose() *
                      metric.value_or(SquareMatrix2::Identity()) * diff)(0, 0));
  }

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lposition Local position (assumed to be in right surface frame)
  /// @param boundaryTolerance boundary check directive
  /// @return boolean indicator for the success of this operation
  virtual bool inside(const Vector2& lposition,
                      const BoundaryTolerance& boundaryTolerance) const {
    using enum BoundaryTolerance::ToleranceMode;

    if (boundaryTolerance.isInfinite()) {
      return true;
    }

    BoundaryTolerance::ToleranceMode mode = boundaryTolerance.toleranceMode();
    bool strictlyInside = inside(lposition);

    if (mode == None) {
      return strictlyInside;
    }

    if (mode == Extend && strictlyInside) {
      return true;
    }

    std::optional<SquareMatrix2> jacobian;
    std::optional<SquareMatrix2> metric;
    if (boundaryTolerance.hasChi2Bound()) {
      SquareMatrix2 j = boundToCartesianJacobian(lposition);
      jacobian = j;
      metric = j.transpose() * boundaryTolerance.asChi2Bound().weight * j;
    } else if (!isCartesian()) {
      jacobian = boundToCartesianJacobian(lposition);
      metric = boundToCartesianMetric(lposition);
    }

    Vector2 closest = closestPoint(lposition, metric);
    Vector2 distance = closest - lposition;

    if (mode == Shrink) {
      return boundaryTolerance.isTolerated(distance, jacobian) &&
             strictlyInside;
    }
    return boundaryTolerance.isTolerated(distance, jacobian);
  }

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param os is the outstream in which the string dump is done
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

}  // namespace Acts
