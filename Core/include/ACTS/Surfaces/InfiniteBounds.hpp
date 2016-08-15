// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// InfiniteBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_INFINITE_BOUNDS_H
#define ACTS_INFINITE_BOUNDS_H 1

#include "ACTS/Surfaces/SurfaceBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

/// @class InfiniteBounds
///
/// templated boundless extension to forward the interface
/// Returns all inside checks to true and can templated for all bounds

class InfiniteBounds : public SurfaceBounds
{
public:
  /// Default Constructor
  InfiniteBounds() {}
  /// Destructor
  ~InfiniteBounds() {}
  /// Return SurfaceBounds type for persistency mainly
  virtual SurfaceBounds::BoundsType
  type() const final
  {
    return SurfaceBounds::Boundless;
  }

  /// Method inside() returns true for any case
  /// ignores input parameters
  /// @return always true
  virtual bool
  inside(const Vector2D&, const BoundaryCheck&) const final
  {
    return true;
  }

  /// Method inside() returns true for loc 0
  /// ignores input parameters
  /// @return always true
  virtual bool
  insideLoc0(const Vector2D& lpos, double tol0 = 0.) const final
  {
    return true;
  }

  /// Method inside() returns true for loc 1
  /// ignores input parameters
  /// @return always true
  virtual bool
  insideLoc1(const Vector2D& lpos, double tol1 = 0.) const final
  {
    return true;
  }

  /// Minimal distance calculation
  /// ignores input parameter
  /// @return always 0. (should be -NaN)
  virtual double
  minDistance(const Vector2D& pos) const final
  {
    return 0.;
  }

  /// Clone method to complete inherited interface
  virtual InfiniteBounds*
  clone() const final
  {
    return new InfiniteBounds();
  }

  /// Output Method for std::ostream
  virtual std::ostream&
  dump(std::ostream& sl) const final;
};

inline std::ostream&
InfiniteBounds::dump(std::ostream& sl) const
{
  sl << "Acts::InfiniteBounds ... boundless surface" << std::endl;
  return sl;
}

static InfiniteBounds s_noBounds = InfiniteBounds();

}  // end of namespace

#endif  // ACTS_INFINITE_BOUNDS_H
