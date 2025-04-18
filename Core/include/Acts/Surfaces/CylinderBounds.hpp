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
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace Acts {

/// @class CylinderBounds
/// @image html CylinderBounds.gif
///
/// Bounds for a cylindrical Surface.
///
/// These bounds may be used for a CylinderSurface
/// In case of bounds for a StraightLineSurface the radius determines the radius
/// within a localPosition
/// is regarded as inside bounds.
///
/// CylinderBounds also enhance the possibility of a cylinder segment with an
/// opening angle @f$ 2\cdot\phi_{half}@f$
/// around an average @f$ \phi @f$ angle @f$ \phi_{ave} @f$.
///
/// CylinderBounds also supports beveled sides defined by an angle.
/// Different angles can be defined on both sides of the cylinder.
/// A positive angle is defined as "extruding" from the defined Zlength,
/// while a negative angle is "intruding" on the Zlength.
/// +    -            -   +
/// ________________________
/// \  |  /          \  |  /
///  \ | /            \ | /
///   \|/______________\|/
///     2 * ZhalfLength
class CylinderBounds : public SurfaceBounds {
 public:
  enum BoundValues : int {
    eR = 0,
    eHalfLengthZ = 1,
    eHalfPhiSector = 2,
    eAveragePhi = 3,
    eBevelMinZ = 4,
    eBevelMaxZ = 5,
    eSize = 6
  };

  CylinderBounds() = delete;

  /// Constructor - full cylinder
  ///
  /// @param r The radius of the cylinder
  /// @param halfZ The half length in z
  /// @param halfPhi The half opening angle
  /// @param avgPhi (optional) The phi value from which the opening angle spans
  /// @param bevelMinZ (optional) The bevel on the negative z side
  /// @param bevelMaxZ (optional) The bevel on the positive z sid The bevel on the positive z side
  CylinderBounds(double r, double halfZ, double halfPhi = M_PI,
                 double avgPhi = 0., double bevelMinZ = 0.,
                 double bevelMaxZ = 0.) noexcept(false)
      : m_values({r, halfZ, halfPhi, avgPhi, bevelMinZ, bevelMaxZ}),
        m_closed(std::abs(halfPhi - M_PI) < s_epsilon) {
    checkConsistency();
  }

  /// Constructor - from fixed size array
  ///
  /// @param values The parameter values
  CylinderBounds(const std::array<double, eSize>& values) noexcept(false)
      : m_values(values),
        m_closed(std::abs(values[eHalfPhiSector] - M_PI) < s_epsilon) {
    checkConsistency();
  }

  ~CylinderBounds() override = default;

  BoundsType type() const final;

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lposition Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  bool inside(const Vector2& lposition,
              const BoundaryCheck& bcheck) const final;

  /// Specialized method for CylinderBounds that checks if a global position
  /// is within the cylinder cover
  ///
  /// @param position is the position in the cylinder frame
  /// @param bcheck is the boundary check directive
  /// @return boolean indicator for operation success
  bool inside3D(const Vector3& position,
                const BoundaryCheck& bcheck = BoundaryCheck(true)) const;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

  /// Returns true for full phi coverage
  bool coversFullAzimuth() const;

  /// Create the bows/circles on either side of the cylinder
  ///
  /// @param trans is the global transform
  /// @param lseg  are the numbero if phi segments
  std::vector<Vector3> createCircles(const Transform3 trans,
                                     std::size_t lseg) const;

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  /// The bound radius, half Z, half phi and average phi
  std::array<double, eSize> m_values;
  /// Indicator if the bounds are closed
  bool m_closed{false};

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);

  /// Helper method to shift into the phi-frame
  /// @param lposition the polar coordinates in the global frame
  Vector2 shifted(const Vector2& lposition) const;

  /// Return the jacobian into the polar coordinate
  ActsMatrix<2, 2> jacobian() const;
};

inline std::vector<double> CylinderBounds::values() const {
  std::vector<double> valvector;
  valvector.insert(valvector.begin(), m_values.begin(), m_values.end());
  return valvector;
}

inline bool CylinderBounds::coversFullAzimuth() const {
  return m_closed;
}

inline void CylinderBounds::checkConsistency() noexcept(false) {
  if (get(eR) <= 0.) {
    throw std::invalid_argument("CylinderBounds: invalid radial setup.");
  }
  if (get(eHalfLengthZ) <= 0.) {
    throw std::invalid_argument("CylinderBounds: invalid length setup.");
  }
  if (get(eHalfPhiSector) <= 0. || get(eHalfPhiSector) > M_PI) {
    throw std::invalid_argument("CylinderBounds: invalid phi sector setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi))) {
    throw std::invalid_argument("CylinderBounds: invalid phi positioning.");
  }
  if (get(eBevelMinZ) != detail::radian_sym(get(eBevelMinZ))) {
    throw std::invalid_argument("CylinderBounds: invalid bevel at min Z.");
  }
  if (get(eBevelMaxZ) != detail::radian_sym(get(eBevelMaxZ))) {
    throw std::invalid_argument("CylinderBounds: invalid bevel at max Z.");
  }
}

}  // namespace Acts
