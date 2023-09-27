// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace Acts {

/// @class DiscTrapezoidBounds
///
/// Class to describe the bounds for a planar DiscSurface.
/// By providing an argument for hphisec, the bounds can
/// be restricted to a phi-range around the center position.

class DiscTrapezoidBounds : public DiscBounds {
 public:
  enum BoundValues : int {
    eHalfLengthXminR = 0,
    eHalfLengthXmaxR = 1,
    eMinR = 2,
    eMaxR = 3,
    eAveragePhi = 4,
    eStereo = 5,
    eSize = 6
  };

  DiscTrapezoidBounds() = delete;

  /// Constructor for a symmetric Trapezoid giving min X length, max X length,
  /// Rmin and R max
  /// @param halfXminR half length in X at min radius
  /// @param halfXmaxR half length in X at maximum radius
  /// @param minR inner radius
  /// @param maxR outer radius
  /// @param avgPhi average phi value
  /// @param stereo optional stero angle applied
  DiscTrapezoidBounds(double halfXminR, double halfXmaxR, double minR,
                      double maxR, double avgPhi = M_PI_2,
                      double stereo = 0.) noexcept(false);

  /// Constructor - from fixed size array
  ///
  /// @param values The parameter values
  DiscTrapezoidBounds(const std::array<double, eSize>& values) noexcept(false)
      : m_values(values) {
    checkConsistency();
  }

  ~DiscTrapezoidBounds() override = default;

  SurfaceBounds::BoundsType type() const final;

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

  ///  This method cheks if the radius given in the LocalPosition is inside
  ///  [rMin,rMax]
  /// if only tol0 is given and additional in the phi sector is tol1 is given
  /// @param lposition is the local position to be checked (in polar
  /// coordinates)
  /// @param bcheck is the boundary check directive
  bool inside(const Vector2& lposition,
              const BoundaryCheck& bcheck = true) const final;

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& sl) const final;

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

  /// This method returns inner radius
  double rMin() const final;

  /// This method returns outer radius
  double rMax() const final;

  /// This method returns the center radius
  double rCenter() const;

  /// This method returns the stereo angle
  double stereo() const;

  /// This method returns the halfPhiSector which is covered by the disc
  double halfPhiSector() const;

  /// This method returns the half length in Y (this is Rmax -Rmin)
  double halfLengthY() const;

  /// Returns true for full phi coverage - obviously false here
  bool coversFullAzimuth() const final;

  /// Checks if this is inside the radial coverage
  /// given the a tolerance
  bool insideRadialBounds(double R, double tolerance = 0.) const final;

  /// Return a reference radius for binning
  double binningValueR() const final;

  /// Return a reference phi for binning
  double binningValuePhi() const final;

  /// This method returns the xy coordinates of the four corners of the
  /// bounds in module coorindates (in xy)
  ///
  /// @param lseg the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note that the number of segments are ignored for this surface
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2> vertices(unsigned int lseg) const final;

 private:
  std::array<double, eSize> m_values;

  /// Dreived maximum y value
  double m_ymax = 0;

  /// Check the input values for consistency, will throw a logic_exception
  /// if consistency is not given
  void checkConsistency() noexcept(false);

  /// Private helper method to convert a local postion
  /// into its Cartesian representation
  ///
  /// @param lposition The local position in polar coordinates
  Vector2 toLocalCartesian(const Vector2& lposition) const;

  /// Jacobian
  /// into its Cartesian representation
  ///
  /// @param lposition The local position in polar coordinates
  ActsMatrix<2, 2> jacobianToLocalCartesian(const Vector2& lposition) const;
};

inline double DiscTrapezoidBounds::rMin() const {
  return get(eMinR);
}

inline double DiscTrapezoidBounds::rMax() const {
  return get(eMaxR);
}

inline double DiscTrapezoidBounds::stereo() const {
  return get(eStereo);
}

inline double DiscTrapezoidBounds::halfPhiSector() const {
  auto minHalfPhi = std::asin(get(eHalfLengthXminR) / get(eMinR));
  auto maxHalfPhi = std::asin(get(eHalfLengthXmaxR) / get(eMaxR));
  return std::max(minHalfPhi, maxHalfPhi);
}

inline double DiscTrapezoidBounds::rCenter() const {
  double rmin = get(eMinR);
  double rmax = get(eMaxR);
  double hxmin = get(eHalfLengthXminR);
  double hxmax = get(eHalfLengthXmaxR);
  auto hmin = std::sqrt(rmin * rmin - hxmin * hxmin);
  auto hmax = std::sqrt(rmax * rmax - hxmax * hxmax);
  return 0.5 * (hmin + hmax);
}

inline double DiscTrapezoidBounds::halfLengthY() const {
  double rmin = get(eMinR);
  double rmax = get(eMaxR);
  double hxmin = get(eHalfLengthXminR);
  double hxmax = get(eHalfLengthXmaxR);
  auto hmin = std::sqrt(rmin * rmin - hxmin * hxmin);
  auto hmax = std::sqrt(rmax * rmax - hxmax * hxmax);
  return 0.5 * (hmax - hmin);
}

inline bool DiscTrapezoidBounds::coversFullAzimuth() const {
  return false;
}

inline bool DiscTrapezoidBounds::insideRadialBounds(double R,
                                                    double tolerance) const {
  return (R + tolerance > get(eMinR) and R - tolerance < get(eMaxR));
}

inline double DiscTrapezoidBounds::binningValueR() const {
  return 0.5 * (get(eMinR) + get(eMaxR));
}

inline double DiscTrapezoidBounds::binningValuePhi() const {
  return get(eAveragePhi);
}

inline std::vector<double> DiscTrapezoidBounds::values() const {
  std::vector<double> valvector;
  valvector.insert(valvector.begin(), m_values.begin(), m_values.end());
  return valvector;
}

inline void DiscTrapezoidBounds::checkConsistency() noexcept(false) {
  if (get(eMinR) < 0. or get(eMaxR) <= 0. or get(eMinR) > get(eMaxR)) {
    throw std::invalid_argument("DiscTrapezoidBounds: invalid radial setup.");
  }
  if (get(eHalfLengthXminR) < 0. or get(eHalfLengthXmaxR) <= 0.) {
    throw std::invalid_argument("DiscTrapezoidBounds: negative length given.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi))) {
    throw std::invalid_argument(
        "DiscTrapezoidBounds: invalid phi positioning.");
  }
}

}  // namespace Acts