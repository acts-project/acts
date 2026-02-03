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
#include "Acts/Surfaces/SurfaceBounds.hpp"

namespace Acts {

/// @class InfiniteBounds
///
/// templated boundless extension to forward the interface
/// Returns all inside checks to true and can templated for all bounds
///
class InfiniteBounds : public SurfaceBounds {
 public:
  /// @copydoc SurfaceBounds::type
  Type type() const final { return Boundless; }

  /// @copydoc SurfaceBounds::isCartesian
  bool isCartesian() const final { return true; }

  /// @copydoc SurfaceBounds::boundToCartesianJacobian
  SquareMatrix2 boundToCartesianJacobian(const Vector2& lposition) const final {
    static_cast<void>(lposition);
    return SquareMatrix2::Identity();
  }

  /// @copydoc SurfaceBounds::boundToCartesianMetric
  SquareMatrix2 boundToCartesianMetric(const Vector2& lposition) const final {
    static_cast<void>(lposition);
    return SquareMatrix2::Identity();
  }

  /// @copydoc SurfaceBounds::values
  std::vector<double> values() const final { return {}; }

  /// @copydoc SurfaceBounds::inside(const Vector2&) const
  bool inside(const Vector2& lposition) const final {
    static_cast<void>(lposition);
    return true;
  }

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const SquareMatrix2& metric) const final {
    static_cast<void>(metric);
    return lposition;
  }

  /// @copydoc SurfaceBounds::inside(const Vector2&, const BoundaryTolerance&) const
  bool inside(const Vector2& lposition,
              const BoundaryTolerance& boundaryTolerance) const final {
    static_cast<void>(lposition);
    static_cast<void>(boundaryTolerance);
    return true;
  }

  /// @copydoc SurfaceBounds::center
  Vector2 center() const final {
    // For infinite bounds, return conceptual center at origin
    return Vector2::Zero();
  }

  /// Output Method for std::ostream
  /// @param os Output stream to write to
  /// @return Reference to the output stream for method chaining
  std::ostream& toStream(std::ostream& os) const final {
    os << "Acts::InfiniteBounds ... boundless surface" << std::endl;
    return os;
  }
};

/// Static instance of InfiniteBounds for convenience
static const InfiniteBounds s_noBounds{};

}  // namespace Acts
