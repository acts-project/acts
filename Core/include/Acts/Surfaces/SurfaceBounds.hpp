// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <ostream>

#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Utilities/Definitions.hpp"

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
class SurfaceBounds {
 public:
  /// @enum BoundsType
  /// This is nested to the SurfaceBounds, as also VolumeBounds will have
  /// Bounds Type.
  enum BoundsType : int {
    eCone = 0,
    eCylinder = 1,
    eDiamond = 2,
    eDisc = 3,
    eEllipse = 5,
    eLine = 6,
    eRectangle = 7,
    eTrapezoid = 8,
    eTriangle = 9,
    eDiscTrapezoid = 10,
    eConvexPolygon = 11,
    eAnnulus = 12,
    eBoundless = 13,
    eOther = 14
  };

  virtual ~SurfaceBounds() = default;

  /// Return the bounds type - for persistency optimization
  ///
  /// @return is a BoundsType enum
  virtual BoundsType type() const = 0;

  /// Access method for bound values, this is a dynamically sized
  /// vector containing the parameters needed to describe these bounds
  ///
  /// @return of the stored values for this SurfaceBounds object
  virtual std::vector<double> values() const = 0;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lposition Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  /// @return boolean indicator for the success of this operation
  virtual bool inside(const Vector2D& lposition,
                      const BoundaryCheck& bcheck) const = 0;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lposition is the local position to check for the distance
  /// @return is a signed distance parameter
  virtual double distanceToBoundary(const Vector2D& lposition) const = 0;

  /// Output Method for std::ostream, to be overloaded by child classes
  ///
  /// @param sl is the outstream in which the string dump is done
  virtual std::ostream& toStream(std::ostream& os) const = 0;
};

inline bool operator==(const SurfaceBounds& lhs, const SurfaceBounds& rhs) {
  if (&lhs == &rhs) {
    return true;
  }
  return (lhs.type() == rhs.type()) && (lhs.values() == rhs.values());
}

inline bool operator!=(const SurfaceBounds& lhs, const SurfaceBounds& rhs) {
  return !(lhs == rhs);
}

inline std::ostream& operator<<(std::ostream& os, const SurfaceBounds& sb) {
  return sb.toStream(os);
}

}  // namespace Acts
