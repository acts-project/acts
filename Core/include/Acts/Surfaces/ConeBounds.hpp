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
#include <cstdlib>
#include <iosfwd>
#include <stdexcept>
#include <vector>

namespace Acts {

///  @class ConeBounds
///
///  Bounds for a conical surface,
///  the opening angle is stored in \f$ \tan(\alpha) \f$ and always positively
/// defined.
///  The cone can open to both sides, steered by \f$ z_min \f$ and \f$ z_max
///  \f$.
///
///  @image html figures/ConeBounds.gif
///

class ConeBounds : public SurfaceBounds {
 public:
  enum BoundValues : int {
    eAlpha = 0,
    eMinZ = 1,
    eMaxZ = 2,
    eHalfPhiSector = 3,
    eAveragePhi = 4,
    eSize = 5
  };

  ConeBounds() = delete;

  /// Constructor - open cone with alpha, by default a full cone
  /// but optionally can make a conical section
  ///
  /// @param alpha is the opening angle of the cone
  /// @param symm is the boolean indicating if the cone is symmetric in +/- z
  /// @param halfphi is the half opening angle (default is pi)
  /// @param avphi is the phi value around which the bounds are opened
  /// (default=0)
  ConeBounds(double alpha, bool symm, double halfphi = M_PI,
             double avphi = 0.) noexcept(false);

  /// Constructor - open cone with alpha, minz and maxz, by
  /// default a full cone but can optionally make it a conical section
  ///
  /// @param alpha is the opening angle of the cone
  /// @param minz cone expanding from minimal z
  /// @param maxz cone expanding to maximal z
  /// @param halfphi is the half opening angle (default is pi)
  /// @param avphi is the phi value around which the bounds are opened
  /// (default=0)
  ConeBounds(double alpha, double minz, double maxz, double halfphi = M_PI,
             double avphi = 0.) noexcept(false);

  /// Constructor - from parameters array
  ///
  /// @param values The parameter array
  ConeBounds(const std::array<double, eSize>& values) noexcept(false);

  ~ConeBounds() override = default;

  BoundsType type() const final;

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  /// inside method for local position
  ///
  /// @param lposition is the local position to be checked
  /// @param bcheck is the boundary check directive
  /// @return is a boolean indicating if the position is inside
  bool inside(const Vector2& lposition,
              const BoundaryCheck& bcheck = true) const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostrea into which the dump is done
  /// @return is the input object
  std::ostream& toStream(std::ostream& sl) const final;

  /// Return the radius at a specific z values
  ///
  /// @param z is the z value for which r is requested
  /// @return is the r value associated with z
  double r(double z) const;

  /// Return tangent of alpha (pre-computed)
  double tanAlpha() const;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

 private:
  std::array<double, eSize> m_values;
  double m_tanAlpha;

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);

  /// Private helper function to shift a local 2D position
  ///
  /// @param lposition The original local position
  Vector2 shifted(const Vector2& lposition) const;
};

inline double ConeBounds::r(double z) const {
  return std::abs(z * m_tanAlpha);
}

inline double ConeBounds::tanAlpha() const {
  return m_tanAlpha;
}

inline std::vector<double> ConeBounds::values() const {
  std::vector<double> valvector;
  valvector.insert(valvector.begin(), m_values.begin(), m_values.end());
  return valvector;
}

inline void ConeBounds::checkConsistency() noexcept(false) {
  if (get(eAlpha) < 0. or get(eAlpha) >= M_PI) {
    throw std::invalid_argument("ConeBounds: invalid open angle.");
  }
  if (get(eMinZ) > get(eMaxZ) or
      std::abs(get(eMinZ) - get(eMaxZ)) < s_epsilon) {
    throw std::invalid_argument("ConeBounds: invalid z range setup.");
  }
  if (get(eHalfPhiSector) < 0. or abs(eHalfPhiSector) > M_PI) {
    throw std::invalid_argument("ConeBounds: invalid phi sector setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi))) {
    throw std::invalid_argument("ConeBounds: invalid phi positioning.");
  }
}

}  // namespace Acts
