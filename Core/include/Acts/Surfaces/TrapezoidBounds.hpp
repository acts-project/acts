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

/// @class TrapezoidBounds
///
/// Bounds for a trapezoidal, planar Surface.
///
/// @image html TrapezoidBounds.gif
///
/// @todo can be speed optimized by calculating kappa/delta and caching it

class TrapezoidBounds : public PlanarBounds {
 public:
  /// @enum BoundValues - for readability
  enum BoundValues {
    bv_minHalfX = 0,
    bv_maxHalfX = 1,
    bv_halfY = 2,
    bv_length = 3
  };

  /// Trapezoid bounds default constructor is forbidden
  TrapezoidBounds() = delete;

  /// Constructor for symmetric Trapezoid
  ///
  /// @param minhalex minimal half length X, definition at negative halflength Y
  /// @param maxhalex maximal half length X, definition at maximum halflength Y
  /// @param haley half length Y - defined at x=0
  TrapezoidBounds(double minhalex, double maxhalex, double haley);

  ~TrapezoidBounds() override;

  TrapezoidBounds* clone() const final;

  BoundsType type() const final;

  std::vector<TDD_real_t> valueStore() const final;

  /// The orientation of the Trapezoid is according to the figure above,
  /// in words: the shorter of the two parallel sides of the trapezoid
  /// intersects
  /// with the negative @f$ y @f$ - axis of the local frame.
  ///
  /// @param lpos is the local position to be checked (Cartesian local frame)
  /// @param bcheck is the boundary check directive
  ///
  /// <br>
  /// The cases are:<br>
  /// (0) @f$ y @f$ or @f$ x @f$ bounds are 0 || 0<br>
  /// (1) the local position is outside @f$ y @f$ bounds <br>
  /// (2) the local position is inside @f$ y @f$ bounds, but outside maximum @f$
  /// x
  /// @f$ bounds  <br>
  /// (3) the local position is inside @f$ y @f$ bounds AND inside minimum @f$ x
  /// @f$ bounds <br>
  /// (4) the local position is inside @f$ y @f$ bounds AND inside maximum @f$ x
  /// @f$ bounds, so that
  /// it depends on the @f$ eta @f$ coordinate
  /// (5) the local position fails test of (4) <br>
  ///
  /// The inside check is done using single equations of straight lines and one
  /// has
  /// to take care if a point
  /// lies on the positive @f$ x @f$ half area(I) or the negative one(II).
  /// Denoting
  /// @f$ |x_{min}| @f$ and
  /// @f$ | x_{max} | @f$ as \c minHalfX respectively \c maxHalfX, such as @f$ |
  /// y_{H} | @f$ as \c halfY,
  /// the equations for the straing lines in (I) and (II) can be written as:<br>
  ///  <br>
  /// - (I):  @f$ y = \kappa_{I} x + \delta_{I} @f$ <br>
  /// - (II): @f$ y = \kappa_{II} x + \delta_{II} @f$ ,<br>
  ///  <br>
  /// where @f$  \kappa_{I} = - \kappa_{II} = 2 \frac{y_{H}}{x_{max} - x_{min}}
  /// @f$
  /// <br>
  /// and   @f$  \delta_{I} = \delta_{II} = - \frac{1}{2}\kappa_{I}(x_{max} +
  /// x_{min}) @f$
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
  /// @note the number of segements is ignored in this representation
  ///
  /// @return vector for vertices in 2D
  std::vector<Vector2D> vertices(unsigned int lseg = 1) const final;

  // Bounding box representation
  const RectangleBounds& boundingBox() const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream& toStream(std::ostream& sl) const final;

  ///  This method returns the minimal halflength in X
  /// (first coordinate of local surface frame)
  double minHalflengthX() const;

  /// This method returns the maximal halflength in X
  /// (first coordinate of local surface frame)
  double maxHalflengthX() const;

  /// This method returns the halflength in Y
  /// (second coordinate of local surface frame)
  double halflengthY() const;

 private:
  double m_minHalfX;
  double m_maxHalfX;
  double m_halfY;
  RectangleBounds m_boundingBox;
};

inline double TrapezoidBounds::minHalflengthX() const {
  return m_minHalfX;
}

inline double TrapezoidBounds::maxHalflengthX() const {
  return m_maxHalfX;
}

inline double TrapezoidBounds::halflengthY() const {
  return m_halfY;
}

}  // namespace Acts