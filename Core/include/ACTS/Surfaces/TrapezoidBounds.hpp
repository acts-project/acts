// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrapezoidBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_TRAPEZOIDBOUNDS_H
#define ACTS_SURFACES_TRAPEZOIDBOUNDS_H 1

#include <cmath>
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/RectangleBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/ParameterDefinitions.hpp"

namespace Acts {

/// @class TrapezoidBounds
///
/// Bounds for a trapezoidal, planar Surface.
///
/// @image html TrapezoidalBounds.gif
///
/// @todo can be speed optimized by calculating kappa/delta and caching it

class TrapezoidBounds : public PlanarBounds
{
public:
  /// @enum BoundValues - for readability
  enum BoundValues {
    bv_minHalfX = 0,
    bv_maxHalfX = 1,
    bv_halfY    = 2,
    bv_length   = 3
  };

  /// Trapezoid bounds default constructor is forbidden
  TrapezoidBounds() = delete;

  /// Constructor for symmetric Trapezoid
  ///
  /// @param minhalex minimal half lenght X, definition at negative halflength Y
  /// @param maxhalex maximal half length X, definition at maximum halflength Y
  /// @param haley half length Y - defined at x=0
  TrapezoidBounds(double minhalex, double maxhalex, double haley);

  /// Copy constructor
  ///
  /// @param trabo are the source bounds for assignment
  TrapezoidBounds(const TrapezoidBounds& trabo)
    : PlanarBounds(trabo), m_boundingBox(0., 0.)
  {
  }

  /// Destructor
  virtual ~TrapezoidBounds();

  /// Virtual constructor
  virtual TrapezoidBounds*
  clone() const final override;

  /// Return the type of the bounds for persistency
  virtual BoundsType
  type() const final
  {
    return SurfaceBounds::Trapezoid;
  }

  /// Assignment operator
  TrapezoidBounds&
  operator=(const TrapezoidBounds& sbo);

  ///  This method returns the minimal halflength in X
  /// (first coordinate of local surface frame)
  double
  minHalflengthX() const;

  /// This method returns the maximal halflength in X
  /// (first coordinate of local surface frame)
  double
  maxHalflengthX() const;

  /// This method returns the halflength in Y
  /// (second coordinate of local surface frame)
  double
  halflengthY() const;

  /// The orientation of the Trapezoid is according to the figure above,
  /// in words: the shorter of the two parallel sides of the trapezoid
  /// intersects
  /// with the negative @f$ y @f$ - axis of the local frame.
  ///
  /// @param lpos is the local position to be checked (carthesian local frame)
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
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  inside(const Vector2D&      lpos,
         const BoundaryCheck& bcheck) const final override;

  /// This method checks inside bounds in loc0
  /// @note loc0/loc1 correspond to the natural coordinates of the surface
  /// @note As loc0/loc1 are correlated the single check doesn't make sense :
  ///       check is done on enclosing Rectangle !
  ///
  /// @param lpos is the local position to be checked
  /// @param tol0 is the tolerance applied
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol0 = 0.) const final override;

  /// This method checks inside bounds in loc0
  /// @note loc0/loc1 correspond to the natural coordinates of the surface
  /// @note As loc0/loc1 are correlated the single check doesn't make sense :
  ///   -> check is done on enclosing Rectangle !
  ///
  /// @param lpos is the local position to be checked
  /// @param tol1 is the tolerance applied
  ///
  /// @return boolean indicator for the success of this operation
  virtual bool
  insideLoc1(const Vector2D& lpos, double tol1 = 0.) const final override;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  ///
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
  /// @param sl is the ostream to be dumped into
  virtual std::ostream&
  dump(std::ostream& sl) const final override;

private:
  RectangleBounds m_boundingBox;  ///< internal bounding box cache
};

inline TrapezoidBounds*
TrapezoidBounds::clone() const
{
  return new TrapezoidBounds(*this);
}

inline double
TrapezoidBounds::minHalflengthX() const
{
  return m_valueStore[TrapezoidBounds::bv_minHalfX];
}

inline double
TrapezoidBounds::maxHalflengthX() const
{
  return m_valueStore[TrapezoidBounds::bv_maxHalfX];
}

inline double
TrapezoidBounds::halflengthY() const
{
  return m_valueStore[TrapezoidBounds::bv_halfY];
}

inline bool
TrapezoidBounds::inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const
{
  return bcheck.isInsidePolygon(lpos, vertices());
}

inline const std::vector<Vector2D>
TrapezoidBounds::vertices() const
{
  // counter-clockwise from bottom-right corner
  return {{minHalflengthX(), -halflengthY()},
          {maxHalflengthX(), halflengthY()},
          {-maxHalflengthX(), halflengthY()},
          {-minHalflengthX(), -halflengthY()}};
}

inline const RectangleBounds&
TrapezoidBounds::boundingBox() const
{
  return m_boundingBox;
}

}  // end of namespace

#endif  // ACTS_SURFACES_TRAPEZOIDBOUNDS_H
