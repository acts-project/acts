// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cmath>
#include <cstdlib>

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class EllipseBounds
///
/// Class to describe the bounds for a planar ellispoid
/// surface.
/// By providing an argument for hphisec, the bounds can
/// be restricted to a phi-range around the center position.
///
/// @image html EllipseBounds.png
///
class EllipseBounds : public PlanarBounds {
 public:
  /// @brief constants for readability
  enum BoundValues {
    bv_rMinX = 0,
    bv_rMinY = 1,
    bv_rMaxX = 2,
    bv_rMaxY = 3,
    bv_averagePhi = 4,
    bv_halfPhiSector = 5,
    bv_length = 6
  };

  /// Deleted default constructor
  EllipseBounds() = delete;

  /// Constructor for full of an ellipsoid disc
  ///
  /// @param minRadius0 is the minimum radius along coordinate 0
  /// @param minRadius1 is the minimum radius along coorindate 1
  /// @param maxRadius0 is the minimum radius at coorindate 0
  /// @param maxRadius1 is the minimum radius at coorindate 1
  /// @param averagePhi average phi (is set to 0. as default)
  /// @param halfPhi    spanning phi sector (is set to pi as default)
  EllipseBounds(double minRadius0, double minRadius1, double maxRadius0,
                double maxRadius1, double averagePhi = 0.,
                double halfPhi = M_PI);

  /// Defaulted destructor
  ~EllipseBounds() override = default;

  /// Clone method for surface cloning
  EllipseBounds* clone() const final;

  /// Type enumeration
  BoundsType type() const final;

  /// Complete value store for persistency
  std::vector<TDD_real_t> valueStore() const final;

  /// This method checks if the point given in the local coordinates is between
  /// two ellipsoids if only tol0 is given and additional in the phi sector is
  /// tol1 is given
  ///
  /// @param lposition Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  bool inside(const Vector2D& lposition,
              const BoundaryCheck& bcheck) const final;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lposition is the local position to check for the distance
  /// @return is a signed distance parameter
  double distanceToBoundary(const Vector2D& lposition) const final;

  /// Return the vertices
  ///
  /// @param lseg the number of segments used to approximate
  /// and eventually curved line, here it refers to the full 2PI Ellipse
  ///
  /// @note the number of segements to may be altered by also providing
  /// the extremas in all direction
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2D> vertices(unsigned int lseg) const final;

  // Bounding box representation
  const RectangleBounds& boundingBox() const final;

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& sl) const final;

  /// This method returns first inner radius
  double rMinX() const;

  /// This method returns second inner radius
  double rMinY() const;

  /// This method returns first outer radius
  double rMaxX() const;

  /// This method returns second outer radius
  double rMaxY() const;

  /// This method returns the average phi
  double averagePhi() const;

  /// This method returns the halfPhiSector which is covered by the disc
  double halfPhiSector() const;

 private:
  double m_rMinX, m_rMinY, m_rMaxX, m_rMaxY, m_avgPhi, m_halfPhi;
  RectangleBounds m_boundingBox;
};

inline double EllipseBounds::rMinX() const {
  return m_rMinX;
}

inline double EllipseBounds::rMinY() const {
  return m_rMinY;
}

inline double EllipseBounds::rMaxX() const {
  return m_rMaxX;
}

inline double EllipseBounds::rMaxY() const {
  return m_rMaxY;
}

inline double EllipseBounds::averagePhi() const {
  return m_avgPhi;
}

inline double EllipseBounds::halfPhiSector() const {
  return m_halfPhi;
}

}  // namespace Acts