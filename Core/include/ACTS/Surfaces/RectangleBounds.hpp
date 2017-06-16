// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RectangleBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_RECTANGLEBOUNDS_H
#define ACTS_SURFACES_RECTANGLEBOUNDS_H 1

#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

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
  /// @param halfX halflength in X
  /// @param halfY halflength in Y
  RectangleBounds(double halfX, double halfY);

  virtual ~RectangleBounds();

  virtual RectangleBounds*
  clone() const final override;

  virtual BoundsType
  type() const final override;

  virtual std::vector<TDD_real_t>
  valueStore() const final override;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  virtual bool
  inside(const Vector2D&      lpos,
         const BoundaryCheck& bcheck) const final override;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  /// @return is a signed distance parameter
  virtual double
  distanceToBoundary(const Vector2D& lpos) const final override;

  /// Return the vertices - or, the points of the extremas
  virtual std::vector<Vector2D>
  vertices() const final override;

  // Bounding box representation
  virtual const RectangleBounds&
  boundingBox() const final override;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream for the dump
  virtual std::ostream&
  dump(std::ostream& sl) const final override;

  /// Return method for the half length in X
  double
  halflengthX() const;

  /// Return method for the half length in Y
  double
  halflengthY() const;

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

}  // end of namespace

#endif  // ACTS_SURFACES_RECTANGLEBOUNDS_H
