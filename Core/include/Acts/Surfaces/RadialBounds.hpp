// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cmath>

#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

namespace Acts {

/// @class RadialBounds
///
/// Class to describe the bounds for a planar DiscSurface.
/// By providing an argument for hphisec, the bounds can
/// be restricted to a phi-range around the center position.
///
/// @image html RadialBounds.gif

class RadialBounds : public DiscBounds {
 public:
  /// enumeration for readability
  enum BoundValues {
    bv_rMin = 0,
    bv_rMax = 1,
    bv_averagePhi = 2,
    bv_halfPhiSector = 3,
    bv_length = 4
  };

  /// Default contructor
  RadialBounds() = delete;

  /// Constructor for full disc of symmetric disc around phi=0
  ///
  /// @param minrad is the inner radius of the disc (0 for full disc)
  /// @param maxrad is the outer radius of the disc
  /// @param hphisec is the half opening angle of the disc (Pi for full angular
  /// coverage)
  RadialBounds(double minrad, double maxrad, double hphisec = M_PI);

  /// Constructor for full disc of symmetric disc around phi!=0
  ///
  /// @param minrad is the inner radius of the disc (0 for full disc)
  /// @param maxrad is the outer radius of the disc
  /// @param avephi is the phi value of the local x-axis in the local 3D frame
  /// @param hphisec is the half opening angle of the disc (Pi for full angular
  /// coverage)
  RadialBounds(double minrad, double maxrad, double avephi, double hphisec);

  /// Defaulted Destructor
  ~RadialBounds() override = default;

  /// Virtual constructor
  RadialBounds* clone() const final;

  /// Return the type enumerator
  SurfaceBounds::BoundsType type() const final;

  /// Value store to be written out
  std::vector<TDD_real_t> valueStore() const final;

  /// For disc surfaces the local position in (r,phi) is checked
  ///
  /// @param lposition local position to be checked
  /// @param bcheck boundary check directive
  ///
  /// @return is a boolean indicating the operation success
  bool inside(const Vector2D& lposition,
              const BoundaryCheck& bcheck) const final;

  /// Minimal distance to boundary calculation
  ///
  /// @param lposition local 2D position in surface coordinate frame
  ///
  /// @return distance to boundary ( > 0 if outside and <=0 if inside)
  double distanceToBoundary(const Vector2D& lposition) const final;

  /// Outstream operator
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream& toStream(std::ostream& sl) const final;

  /// Return method for inner Radius
  double rMin() const;

  /// Return method for outer Radius
  double rMax() const;

  /// Return method for the central phi value
  ///(i.e. phi value of x-axis of local 3D frame)
  double averagePhi() const;

  /// Return method for the half phi span
  double halfPhiSector() const;

  /// Returns true for full phi coverage
  bool coversFullAzimuth() const final;

  /// Checks if this is inside the radial coverage
  /// given the a tolerance
  bool insideRadialBounds(double R, double tolerance = 0.) const final;

  /// Return a reference radius for binning
  double binningValueR() const final;

  /// Return a reference radius for binning
  double binningValuePhi() const final;

 private:
  double m_rMin, m_rMax, m_avgPhi, m_halfPhi;

  /// Private helper method to shift a local position
  /// within the bounds
  ///
  /// @param lposition The local position in polar coordinates
  Vector2D shifted(const Vector2D& lposition) const;

  /// This method returns the xy coordinates of vertices along
  /// the radial bounds
  ///
  /// @param lseg the number of segments used to approximate
  /// and eventually curved line
  ///
  /// @note that the extremas are given, which may slightly alter the
  /// number of segments returned
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2D> vertices(unsigned int lseg) const;
};

inline double RadialBounds::rMin() const {
  return m_rMin;
}

inline double RadialBounds::rMax() const {
  return m_rMax;
}

inline double RadialBounds::averagePhi() const {
  return m_avgPhi;
}

inline double RadialBounds::halfPhiSector() const {
  return m_halfPhi;
}

inline bool RadialBounds::coversFullAzimuth() const {
  return (m_halfPhi == M_PI);
}

inline bool RadialBounds::insideRadialBounds(double R, double tolerance) const {
  return (R + tolerance > m_rMin and R - tolerance < m_rMax);
}

inline double RadialBounds::binningValueR() const {
  return 0.5 * (m_rMin + m_rMax);
}

inline double RadialBounds::binningValuePhi() const {
  return m_avgPhi;
}

}  // namespace Acts