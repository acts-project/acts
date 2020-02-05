// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscTrapezoidalBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <cmath>

#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

///
/// @class DiscTrapezoidalBounds
///
/// Class to describe the bounds for a planar DiscSurface.
/// By providing an argument for hphisec, the bounds can
/// be restricted to a phi-range around the center position.
///

class DiscTrapezoidalBounds : public DiscBounds {
 public:
  /// @enum BoundValues
  /// enumeration for readability
  enum BoundValues {
    bv_rMin = 0,
    bv_rMax = 1,
    bv_minHalfX = 2,
    bv_maxHalfX = 3,
    bv_averagePhi = 4,
    bv_stereo = 5,
    bv_length = 6
  };

  DiscTrapezoidalBounds() = delete;

  /// Constructor for a symmetric Trapezoid giving min X length, max X length,
  /// Rmin and R max
  /// @param minhalfx half length in X at min radius
  /// @param maxhalfx half length in X at maximum radius
  /// @param maxR outer radius
  /// @param minR inner radius
  /// @param avephi average phi value
  /// @param stereo optional stero angle applied
  DiscTrapezoidalBounds(double minhalfx, double maxhalfx, double maxR,
                        double minR, double avephi = M_PI_2,
                        double stereo = 0.);

  ~DiscTrapezoidalBounds() override;

  DiscTrapezoidalBounds* clone() const final;

  SurfaceBounds::BoundsType type() const final;

  std::vector<TDD_real_t> valueStore() const final;

  ///  This method cheks if the radius given in the LocalPosition is inside
  ///  [rMin,rMax]
  /// if only tol0 is given and additional in the phi sector is tol1 is given
  /// @param lposition is the local position to be checked (in polar
  /// coordinates)
  /// @param bcheck is the boundary check directive
  bool inside(const Vector2D& lposition,
              const BoundaryCheck& bcheck = true) const final;

  /// Minimal distance to boundary
  /// @param lposition is the local position to be checked (in polar
  /// coordinates)
  /// @return is the minimal distance ( > 0 if outside and <=0 if inside)
  double distanceToBoundary(const Vector2D& lposition) const final;

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& sl) const final;

  /// This method returns inner radius
  double rMin() const;

  /// This method returns outer radius
  double rMax() const;

  /// This method returns the average phi
  double averagePhi() const;

  /// This method returns the center radius
  double rCenter() const;

  /// This method returns the stereo angle
  double stereo() const;

  /// This method returns the halfPhiSector which is covered by the disc
  double halfPhiSector() const;

  /// This method returns the minimal halflength in X
  double minHalflengthX() const;

  /// This method returns the maximal halflength in X
  double maxHalflengthX() const;

  /// This method returns the halflength in Y (this is Rmax -Rmin)
  double halflengthY() const;

  /// Returns true for full phi coverage - obviously false here
  bool coversFullAzimuth() const final;

  /// Checks if this is inside the radial coverage
  /// given the a tolerance
  bool insideRadialBounds(double R, double tolerance = 0.) const final;

  /// Return a reference radius for binning
  double binningValueR() const final;

  /// Return a reference phi for binning
  double binningValuePhi() const final;

 private:
  double m_rMin, m_rMax, m_minHalfX, m_maxHalfX, m_avgPhi;
  double m_stereo;  // TODO 2017-04-09 msmk: what is this good for?

  /// Private helper method to convert a local postion
  /// into its Cartesian representation
  ///
  /// @param lposition The local position in polar coordinates
  Vector2D toLocalCartesian(const Vector2D& lposition) const;

  /// Jacobian
  /// into its Cartesian representation
  ///
  /// @param lposition The local position in polar coordinates
  ActsMatrixD<2, 2> jacobianToLocalCartesian(const Vector2D& lposition) const;
};

inline double DiscTrapezoidalBounds::rMin() const {
  return m_rMin;
}

inline double DiscTrapezoidalBounds::rMax() const {
  return m_rMax;
}

inline double DiscTrapezoidalBounds::minHalflengthX() const {
  return m_minHalfX;
}

inline double DiscTrapezoidalBounds::maxHalflengthX() const {
  return m_maxHalfX;
}

inline double DiscTrapezoidalBounds::averagePhi() const {
  return m_avgPhi;
}

inline double DiscTrapezoidalBounds::stereo() const {
  return m_stereo;
}

inline double DiscTrapezoidalBounds::halfPhiSector() const {
  auto minHalfPhi = std::asin(m_minHalfX / m_rMin);
  auto maxHalfPhi = std::asin(m_maxHalfX / m_rMax);
  return std::max(minHalfPhi, maxHalfPhi);
}

inline double DiscTrapezoidalBounds::rCenter() const {
  auto hmin = std::sqrt(m_rMin * m_rMin - m_minHalfX * m_minHalfX);
  auto hmax = std::sqrt(m_rMax * m_rMax - m_maxHalfX * m_maxHalfX);
  return (hmin + hmax) / 2.0;
}

inline double DiscTrapezoidalBounds::halflengthY() const {
  auto hmin = std::sqrt(m_rMin * m_rMin - m_minHalfX * m_minHalfX);
  auto hmax = std::sqrt(m_rMax * m_rMax - m_maxHalfX * m_maxHalfX);
  return (hmax - hmin) / 2.0;
}

inline bool DiscTrapezoidalBounds::coversFullAzimuth() const {
  return false;
}

inline bool DiscTrapezoidalBounds::insideRadialBounds(double R,
                                                      double tolerance) const {
  return (R + tolerance > m_rMin and R - tolerance < m_rMax);
}

inline double DiscTrapezoidalBounds::binningValueR() const {
  return 0.5 * (m_rMin + m_rMax);
}

inline double DiscTrapezoidalBounds::binningValuePhi() const {
  return m_avgPhi;
}

}  // namespace Acts