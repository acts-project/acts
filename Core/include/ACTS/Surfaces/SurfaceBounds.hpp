// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_SURFACEBOUNDS_H
#define ACTS_SURFACES_SURFACEBOUNDS_H 1

// STD/STL
#include <iomanip>
#include <iostream>

#include "ACTS/Surfaces/BoundaryCheck.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

/**
 @class SurfaceBounds

 Abstract base class for surface bounds to be specified.

 Surface bounds provide:
 - inside() checks
 - the BoundsType return type to avoid dynamic casting
 - an initCache() method

 */

class SurfaceBounds
{
public:
  /** @enum BoundsType

      This enumerator simplifies the persistency,
      by saving a dynamic_cast to happen.

      Other is reserved for the GeometrySurfaces implementation.

    */
  enum BoundsType {
    Cone             = 0,
    Cylinder         = 1,
    Diamond          = 2,
    Disc             = 3,
    Ellipse          = 5,
    Rectangle        = 6,
    RotatedTrapezoid = 7,
    Trapezoid        = 8,
    Triangle         = 9,
    DiscTrapezoidal  = 10,
    Other            = 11
  };

  /**Default Constructor*/
  SurfaceBounds() {}
  /**Destructor*/
  virtual ~SurfaceBounds() {}
  /** clone() method to make deep copy in Surface copy constructor and for
    assigment operator
    of the Surface class.*/
  virtual SurfaceBounds*
  clone() const = 0;

  /**Equality operator*/
  virtual bool
  operator==(const SurfaceBounds& sb) const = 0;

  /**Non-Equality operator*/
  bool
  operator!=(const SurfaceBounds& sb) const;

  /** Return the bounds type - for persistency optimization */
  virtual BoundsType
  type() const = 0;

  /** Each Bounds has a method inside, which checks if a LocalPosition is inside
     the bounds.
      Inside can be called without/with tolerances. */
  virtual bool
  inside(const Vector2D& locpo, double tol1 = 0., double tol2 = 0.) const = 0;
  virtual bool
  inside(const Vector2D& locpo, const BoundaryCheck& bchk) const = 0;

  /** Extend the interface to for single inside Loc 1 / Loc2 tests
     - loc1/loc2 correspond to the natural coordinates of the surface */
  virtual bool
  insideLoc1(const Vector2D& locpo, double tol1 = 0.) const = 0;

  /** Extend the interface to for single inside Loc 1 / Loc2 tests
     - loc1/loc2 correspond to the natural coordinates of the surface */
  virtual bool
  insideLoc2(const Vector2D& locpo, double tol2 = 0.) const = 0;

  /** Minimal distance to boundary ( > 0 if outside and <=0 if inside) */
  virtual double
  minDistance(const Vector2D& pos) const = 0;

  /** Interface method for the maximal extension or the radius*/
  virtual double
  r() const = 0;

  /** Output Method for std::ostream, to be overloaded by child classes */
  virtual std::ostream&
  dump(std::ostream& sl) const = 0;

protected:
  /** Swap method to be called from DiscBounds or TrapezoidalBounds */
  void
  swap(double& b1, double& b2);

  /** virtual initCache method for object persistency */
  virtual void
  initCache()
  {
  }
};

inline void
SurfaceBounds::swap(double& b1, double& b2)
{
  double tmp = b1;
  b1         = b2;
  b2         = tmp;
}

inline bool
SurfaceBounds::operator!=(const SurfaceBounds& sb) const
{
  return !((*this) == sb);
}

/**Overload of << operator for std::ostream for debug output*/
std::ostream&
operator<<(std::ostream& sl, const SurfaceBounds& sb);

}  // end of namespace

#endif  // ACTS_SURFACES_SURFACEBOUNDS_H
