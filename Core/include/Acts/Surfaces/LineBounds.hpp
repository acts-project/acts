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

#include <array>
#include <iosfwd>
#include <vector>

namespace Acts {

/// @class LineBounds
///
/// Bounds for a LineSurface.
class LineBounds : public SurfaceBounds {
 public:
  /// @enum BoundValues
  /// Enumeration for the bound values
  enum BoundValues : int { eR = 0, eHalfLengthZ = 1, eSize = 2 };

  /// Constructor
  ///
  /// @param r The radius of the line
  /// @param halfZ The half length in z
  explicit LineBounds(double r, double halfZ) noexcept(false)
      : m_values({r, halfZ}) {
    checkConsistency();
  }

  /// Constructor - from fixed size array
  ///
  /// @param values The bound values stored in a fixed size array
  explicit LineBounds(const std::array<double, eSize>& values) noexcept(false)
      : m_values(values) {
    checkConsistency();
  }

  /// @copydoc SurfaceBounds::type
  BoundsType type() const final { return Line; }

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

  /// Return the bound values as dynamically sized vector
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// @copydoc SurfaceBounds::inside
  bool inside(const Vector2& lposition) const final;

  /// @copydoc SurfaceBounds::closestPoint
  Vector2 closestPoint(const Vector2& lposition,
                       const SquareMatrix2& metric) const final;

  using SurfaceBounds::inside;

  /// @copydoc SurfaceBounds::center
  /// @note For LineBounds: returns (0,0) since bounds are symmetric around origin
  Vector2 center() const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream to be dumped into
  /// @return Reference to the output stream for method chaining
  std::ostream& toStream(std::ostream& sl) const final;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  /// @return The bound value for the specified parameter
  double get(BoundValues bValue) const { return m_values[bValue]; }

 private:
  std::array<double, eSize> m_values;

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);
};

}  // namespace Acts
