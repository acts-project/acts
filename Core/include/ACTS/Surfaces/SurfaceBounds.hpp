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

/// @class SurfaceBounds
///
/// Abstract base class for surface bounds to be specified.
///
/// Surface bounds provide:
/// - inside() checks
/// - the BoundsType return type to avoid dynamic casting
///
/// @todo for easy persistency access, force Constructor from ValueStore
///
class SurfaceBounds
{
public:
  /// @enum BoundsType
  ///
  /// This enumerator simplifies the persistency,
  /// by saving a dynamic_cast to happen.
  ///
  enum BoundsType {
    Cone             = 0,
    Cylinder         = 1,
    Diamond          = 2,
    Disc             = 3,
    Ellipse          = 5,
    Line             = 6,
    Rectangle        = 7,
    RotatedTrapezoid = 8,
    Trapezoid        = 9,
    Triangle         = 10,
    DiscTrapezoidal  = 11,
    Boundless        = 12,
    Other            = 12
  };

  /// Default Constructor
  /// @param sSize is the size of the data store
  /// @note the value Store is initialized to the given size
  SurfaceBounds(size_t sSize = 0) : m_valueStore(sSize, 0.) {}
  /// Copy constructor
  /// It copies the value store
  SurfaceBounds(const SurfaceBounds& sb) : m_valueStore(sb.m_valueStore) {}
  /// Destructor
  virtual ~SurfaceBounds() {}
  /// clone() method to make deep copy in Surface copy constructor and for
  /// assigment operator of the Surface class
  virtual SurfaceBounds*
  clone() const = 0;

  /// Assignment operator
  SurfaceBounds&
  operator=(const SurfaceBounds& sb);

  /// Comparison (equality) operator
  /// checks first on the pointer equality
  /// then it cheks on the type
  /// lastly it checks on the data store
  virtual bool
  operator==(const SurfaceBounds& sb) const;

  /// Comparison (non-equality) operator
  /// checks first on the pointer equality, inverts operator==
  bool
  operator!=(const SurfaceBounds& sb) const;

  /// Return the bounds type - for persistency optimization
  virtual BoundsType
  type() const = 0;

  /// Access method for bound variable store
  /// @return of the stored values for the boundary object
  virtual const std::vector<TDD_real_t>&
  valueStore() const;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bchk boundary check directive
  /// @return boolean indicator for the success of this operation
  virtual bool
  inside(const Vector2D& lpos, const BoundaryCheck& bchk) const = 0;

  /// Inside check for the bounds object with tolerance
  /// checks for first coordinate only.
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol1 tolerance parameter
  /// @return boolean indicator for the success of this operation
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol1 = 0.) const = 0;

  /// Inside check for the bounds object with tolerance
  /// checks for second coordinate only.
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param tol2 tolerance parameter
  /// @return boolean indicator for the success of this operation
  virtual bool
  insideLoc1(const Vector2D& lpos, double tol2 = 0.) const = 0;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  virtual double
  minDistance(const Vector2D& pos) const = 0;

  /// Output Method for std::ostream, to be overloaded by child classes
  virtual std::ostream&
  dump(std::ostream& sl) const = 0;

protected:
  std::vector<TDD_real_t> m_valueStore;  ///< internal data store
};

inline bool
SurfaceBounds::operator!=(const SurfaceBounds& sb) const
{
  return !((*this) == sb);
}

inline const std::vector<TDD_real_t>&
SurfaceBounds::valueStore() const
{
  return m_valueStore;
}

std::ostream&
operator<<(std::ostream& sl, const SurfaceBounds& sb);

}  // end of namespace

#endif  // ACTS_SURFACES_SURFACEBOUNDS_H
