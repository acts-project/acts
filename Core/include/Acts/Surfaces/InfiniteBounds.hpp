// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

namespace Acts {

/// @class InfiniteBounds
///
/// templated boundless extension to forward the interface
/// Returns all inside checks to true and can templated for all bounds

class InfiniteBounds : public SurfaceBounds {
 public:
  SurfaceBounds::BoundsType type() const final { return eBoundless; }

  bool isCartesian() const final { return true; }

  SquareMatrix2 boundToCartesianJacobian(const Vector2& lposition) const final {
    (void)lposition;
    return SquareMatrix2::Identity();
  }

  SquareMatrix2 cartesianToBoundJacobian(const Vector2& lposition) const final {
    (void)lposition;
    return SquareMatrix2::Identity();
  }

  SquareMatrix2 boundToCartesianMetric(const Vector2& lposition) const final {
    (void)lposition;
    return SquareMatrix2::Identity();
  }

  std::vector<double> values() const final { return {}; }

  bool inside(const Vector2& /*lposition*/) const final { return true; }

  Vector2 closestPoint(
      const Vector2& lposition,
      const std::optional<SquareMatrix2>& /*metric*/) const final {
    return lposition;
  }

  /// Method inside() returns true for any case
  ///
  /// ignores input parameters
  ///
  /// @return always true
  bool inside(const Vector2& /*lposition*/,
              const BoundaryTolerance& /*boundaryTolerance*/) const final {
    return true;
  }

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& os) const final {
    os << "Acts::InfiniteBounds ... boundless surface" << std::endl;
    return os;
  }
};

static const InfiniteBounds s_noBounds{};

}  // namespace Acts
