// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cmath>

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {

///
/// @class DiamondBounds
///
/// Bounds for a double trapezoidal ("diamond"), planar Surface.
///
class DiamondBounds : public PlanarBounds {
 public:
  /// @enum BoundValues for better readability
  enum BoundValues {
    bv_x1 = 0,
    bv_x2 = 1,
    bv_x3 = 2,
    bv_y1 = 3,
    bv_y2 = 4,
    bv_length = 5
  };

  /// Constructor for convex hexagon symmetric about the y axis
  ///
  /// @param x1 is the halflength in x at minimal y
  /// @param x2 is the halflength in x at y = 0
  /// @param x3 is the halflength in x at maximal y
  /// @param y1 is the halflength into y < 0
  /// @param y2 is the halflength into y > 0
  ///
  /// @image html DiamondBounds.svg
  DiamondBounds(double x1, double x2, double x3, double y1, double y2);

  /// Defaulted desctructor
  ~DiamondBounds() override = default;

  /// Virtual constructor
  DiamondBounds* clone() const final;

  /// Enumeration type
  BoundsType type() const final;

  /// The value store for persistency
  std::vector<TDD_real_t> valueStore() const final;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
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
  /// and eventually curved line
  ///
  /// @note the number of segements is ignored for this representation
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2D> vertices(unsigned int lseg = 1) const final;

  // Bounding box representation
  const RectangleBounds& boundingBox() const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream in which it is dumped
  std::ostream& toStream(std::ostream& sl) const final;

  /// This method returns the halflength in X at minimal Y
  /// (first coordinate of local surface frame)
  double x1() const;

  /// This method returns the (maximal) halflength in X
  /// (first coordinate of local surface frame)
  double x2() const;

  /// This method returns the halflength in X at maximal Y
  /// (first coordinate of local surface frame)
  double x3() const;

  /// This method returns the halflength in Y of trapezoid at negative Y
  double y1() const;

  /// This method returns the halflength in Y of trapezoid at positive Y
  double y2() const;

 private:
  double m_x1, m_x2, m_x3;
  double m_y1, m_y2;
  RectangleBounds m_boundingBox;  ///< internal bounding box cache
};

inline double DiamondBounds::x1() const {
  return m_x1;
}

inline double DiamondBounds::x2() const {
  return m_x2;
}

inline double DiamondBounds::x3() const {
  return m_x3;
}

inline double DiamondBounds::y1() const {
  return m_y1;
}

inline double DiamondBounds::y2() const {
  return m_y2;
}

}  // namespace Acts