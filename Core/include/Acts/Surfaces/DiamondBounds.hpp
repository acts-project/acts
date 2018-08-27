// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiamondBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <cmath>

#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

///
/// @class DiamondBounds
///
/// Bounds for a double trapezoidal ("diamond"), planar Surface.
///
class DiamondBounds : public PlanarBounds
{
public:
  /// @enum BoundValues for better readability
  enum BoundValues {
    bv_minHalfX = 0,
    bv_medHalfX = 1,
    bv_maxHalfX = 2,
    bv_halfY1   = 3,
    bv_halfY2   = 4,
    bv_length   = 5
  };

  /// Constructor for convex hexagon symmetric about the y axis
  ///
  /// @param minhalex is the halflength in x at minimal y
  /// @param medhalex is the halflength in x at y = 0
  /// @param maxhalex is the halflength in x at maximal y
  /// @param haley1 is the halflength into y < 0
  /// @param haley2 is the halflength into y > 0
  DiamondBounds(double minhalex,
                double medhalex,
                double maxhalex,
                double haley1,
                double haley2);

  /// Constructor which accepts @c variant_data
  ///
  /// @param vardata the @c variant_data to build from
  DiamondBounds(const variant_data& vardata);

  ~DiamondBounds() override;

  DiamondBounds*
  clone() const final;

  BoundsType
  type() const final;

  std::vector<TDD_real_t>
  valueStore() const final;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  bool
  inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const final;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  /// @return is a signed distance parameter
  double
  distanceToBoundary(const Vector2D& lpos) const final;

  /// Return the vertices - or, the points of the extremas
  std::vector<Vector2D>
  vertices() const final;

  // Bounding box representation
  const RectangleBounds&
  boundingBox() const final;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream in which it is dumped
  std::ostream&
  dump(std::ostream& sl) const final;

  /// This method returns the halflength in X at minimal Y
  /// (first coordinate of local surface frame)
  double
  minHalflengthX() const;

  /// This method returns the (maximal) halflength in X
  /// (first coordinate of local surface frame)
  double
  medHalflengthX() const;

  /// This method returns the halflength in X at maximal Y
  /// (first coordinate of local surface frame)
  double
  maxHalflengthX() const;

  /// This method returns the halflength in Y of trapezoid at negative Y
  double
  halflengthY1() const;

  /// This method returns the halflength in Y of trapezoid at positive Y
  double
  halflengthY2() const;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const override;

private:
  double          m_minHalfX, m_medHalfX, m_maxHalfX;
  double          m_minY, m_maxY;
  RectangleBounds m_boundingBox;  ///< internal bounding box cache
};

inline double
DiamondBounds::minHalflengthX() const
{
  return m_minHalfX;
}

inline double
DiamondBounds::medHalflengthX() const
{
  return m_medHalfX;
}

inline double
DiamondBounds::maxHalflengthX() const
{
  return m_maxHalfX;
}

inline double
DiamondBounds::halflengthY1() const
{
  return m_minY;
}

inline double
DiamondBounds::halflengthY2() const
{
  return m_maxY;
}

}  // namespace