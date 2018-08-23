// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <ostream>

#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

/// @class SurfaceBounds
///
/// Interface for surface bounds.
///
/// Surface bounds provide:
/// - inside() checks
/// - distance to boundary calculations
/// - the BoundsType and a set of parameters to simplify persistency
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

  virtual ~SurfaceBounds() = default;

  /// clone() method to make deep copy in Surface copy constructor and for
  /// assigment operator of the Surface class
  ///
  /// @return is a newly created object
  virtual SurfaceBounds*
  clone() const = 0;

  /// Return the bounds type - for persistency optimization
  ///
  /// @return is a BoundsType enum
  virtual BoundsType
  type() const = 0;

  /// Access method for bound variable store
  ///
  /// @return of the stored values for the boundary object
  virtual std::vector<TDD_real_t>
  valueStore() const = 0;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lpos Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  virtual bool
  inside(const Vector2D& lpos, const BoundaryCheck& bcheck) const = 0;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lpos is the local position to check for the distance
  /// @return is a signed distance parameter
  virtual double
  distanceToBoundary(const Vector2D& lpos) const = 0;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the outstream in which the string dump is done
  virtual std::ostream&
  dump(std::ostream& os) const = 0;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  virtual variant_data
  toVariantData() const = 0;
};

inline bool
operator==(const SurfaceBounds& lhs, const SurfaceBounds& rhs)
{
  if (&lhs == &rhs) {
    return true;
  }
  return (lhs.type() == rhs.type()) && (lhs.valueStore() == rhs.valueStore());
}

inline bool
operator!=(const SurfaceBounds& lhs, const SurfaceBounds& rhs)
{
  return !(lhs == rhs);
}

inline std::ostream&
operator<<(std::ostream& os, const SurfaceBounds& sb)
{
  return sb.dump(os);
}

}  // namespace Acts
