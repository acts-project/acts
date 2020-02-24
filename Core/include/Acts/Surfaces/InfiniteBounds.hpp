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

/// @class InfiniteBounds
///
/// templated boundless extension to forward the interface
/// Returns all inside checks to true and can templated for all bounds

class InfiniteBounds : public SurfaceBounds {
 public:
  InfiniteBounds() = default;
  ~InfiniteBounds() override = default;

  InfiniteBounds* clone() const final { return new InfiniteBounds(); }

  SurfaceBounds::BoundsType type() const final {
    return SurfaceBounds::Boundless;
  }

  std::vector<TDD_real_t> valueStore() const final { return {}; }

  /// Method inside() returns true for any case
  ///
  /// ignores input parameters
  ///
  /// @return always true
  bool inside(const Vector2D& /*lposition*/,
              const BoundaryCheck& /*bcheck*/) const final {
    return true;
  }

  /// Minimal distance calculation
  /// ignores input parameter
  /// @return always 0. (should be -NaN)
  double distanceToBoundary(const Vector2D& /*position*/) const final {
    return 0;
  }

  /// Output Method for std::ostream
  std::ostream& toStream(std::ostream& os) const final {
    os << "Acts::InfiniteBounds ... boundless surface" << std::endl;
    return os;
  }
};

static const InfiniteBounds s_noBounds{};

}  // namespace Acts
