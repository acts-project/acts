// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @class LineBounds
///
/// Bounds for a LineSurface.
///

class LineBounds : public SurfaceBounds {
 public:
  /// @enum BoundValues for readablility
  /// nested enumeration object
  enum BoundValues { bv_radius = 0, bv_halfZ = 1, bv_length = 2 };

  /// Constructor
  ///
  /// @param radius is the radius of the cylinder, default = 0.
  /// @param halez is the half length in z, defualt = 0.
  LineBounds(double radius = 0., double halez = 0.);

  /// Defaulted destructor
  ~LineBounds() override = default;

  /// Virtual constructor
  LineBounds* clone() const final;

  /// Type enumeration
  BoundsType type() const final;

  /// Return the intrinsic values
  std::vector<TDD_real_t> valueStore() const final;

  /// Inside check for the bounds object driven by the boundary check directive
  /// Each Bounds has a method inside, which checks if a LocalPosition is inside
  /// the bounds  Inside can be called without/with tolerances.
  ///
  /// @param lposition Local position (assumed to be in right surface frame)
  /// @param bcheck boundary check directive
  ///
  /// @return boolean indicator for the success of this operation
  bool inside(const Vector2D& lposition,
              const BoundaryCheck& bcheck) const final;

  /// Minimal distance to boundary ( > 0 if outside and <=0 if inside)
  ///
  /// @param lposition is the local position to check for the distance
  ///
  /// @return is a signed distance parameter
  double distanceToBoundary(const Vector2D& lposition) const final;

  /// This method returns the radius
  virtual double r() const;

  /// This method returns the halflengthZ
  double halflengthZ() const;

  /// Output Method for std::ostream
  ///
  /// @param sl is the ostream to be dumped into
  std::ostream& toStream(std::ostream& sl) const final;

 private:
  double m_radius, m_halfZ;
};

inline double LineBounds::r() const {
  return m_radius;
}

inline double LineBounds::halflengthZ() const {
  return m_halfZ;
}

}  // namespace Acts