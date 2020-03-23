// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace Acts {

/// @class EllipseBounds
///
/// Class to describe the bounds for a planar ellispoid
/// surface.
///
/// By providing an argument for hphisec, the bounds can
/// be restricted to a phi-range around the center position.
class EllipseBounds : public PlanarBounds {
 public:
  enum BoundValues {
    eMinR0 = 0,
    eMinR1 = 1,
    eMaxR0 = 2,
    eMaxR1 = 3,
    eAveragePhi = 4,
    eHalfPhiSector = 5,
    eSize = 6
  };

  EllipseBounds() = delete;

  /// Constructor for full of an ellipsoid disc
  ///
  /// @param minR0 is the minimum radius along coordinate 0
  /// @param minR1 is the minimum radius along coorindate 1
  /// @param maxR0 is the minimum radius at coorindate 0
  /// @param maxR1 is the minimum radius at coorindate 1
  /// @param averagePhi average phi (is set to 0. as default)
  /// @param halfPhi spanning phi sector (is set to pi as default)
  EllipseBounds(double minR0, double minR1, double maxR0, double maxR1,
                double averagePhi = 0., double halfPhi = M_PI);

  /// Constructor - from fixed size array
  ///
  /// @param values The parameter values
  EllipseBounds(const std::array<double, eSize>& values)
      : m_values(values), m_boundingBox(values[eMaxR0], values[eMaxR1]) {}

  ~EllipseBounds() override = default;

  EllipseBounds* clone() const final;

  BoundsType type() const final;

  /// Return the bound values as dynamically sized vector
  ///
  /// @return this returns a copy of the internal values
  std::vector<double> values() const final;

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

  /// Access to the bound values
  /// @param bValue the class nested enum for the array access
  double get(BoundValues bValue) const { return m_values[bValue]; }

 private:
  std::array<double, eSize> m_values;
  RectangleBounds m_boundingBox;
};

inline std::vector<double> EllipseBounds::values() const {
  std::vector<double> valvector;
  valvector.insert(valvector.begin(), m_values.begin(), m_values.end());
  return valvector;
}

}  // namespace Acts