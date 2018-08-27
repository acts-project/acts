// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RectangleBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

/// @class RectangleBounds
///
/// Bounds for a rectangular, planar surface.
/// The two local coordinates Acts::eLOC_X, Acts::eLOC_Y are for legacy reasons
/// also called @f$ phi @f$ respectively @f$ eta @f$. The orientation
/// with respect to the local surface framce can be seen in the attached
/// illustration.
///
/// @image html RectangularBounds.gif

class RectangleBounds : public PlanarBounds
{
public:
  /// @enum BoundValues for readability
  enum BoundValues { bv_halfX = 0, bv_halfY = 1, bv_length = 2 };

  RectangleBounds() = delete;

  /// Constructor with halflength in x and y
  ///
  /// @param halex halflength in X
  /// @param haley halflength in Y
  RectangleBounds(double halex, double haley);

  /// Constructor which accepts @c variant_data
  ///
  /// @param vardata the @c variant_data to build from
  RectangleBounds(const variant_data& vardata);

  ~RectangleBounds() override;

  RectangleBounds*
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
  /// @param sl is the ostream for the dump
  std::ostream&
  dump(std::ostream& sl) const final;

  /// Return method for the half length in X
  double
  halflengthX() const;

  /// Return method for the half length in Y
  double
  halflengthY() const;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const override;

private:
  double m_halfX, m_halfY;
};

inline double
RectangleBounds::halflengthX() const
{
  return m_halfX;
}

inline double
RectangleBounds::halflengthY() const
{
  return m_halfY;
}

inline SurfaceBounds::BoundsType
RectangleBounds::type() const
{
  return SurfaceBounds::Rectangle;
}

}  // namespace