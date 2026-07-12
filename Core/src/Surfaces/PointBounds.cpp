// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PointBounds.hpp"

#include "Acts/Utilities/detail/OstreamStateGuard.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace Acts {

std::vector<double> PointBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

void PointBounds::checkConsistency() noexcept(false) {
  if (get(eR) < 0.) {
    throw std::invalid_argument("PointBounds: negative radius.");
  }
}

bool PointBounds::inside(const Vector2& lposition) const {
  double r = get(eR);
  return lposition.squaredNorm() <= r * r;
}

Vector2 PointBounds::closestPoint(const Vector2& lposition,
                                  const SquareMatrix2& metric) const {
  static_cast<void>(metric);
  // The bounds are a disc of radius R; the closest point on the boundary
  // circle is the radial projection of the local position. For a position at
  // the origin any boundary point is equidistant, pick (R, 0).
  double r = get(eR);
  double norm = lposition.norm();
  if (norm == 0.) {
    return Vector2(r, 0.);
  }
  return lposition * (r / norm);
}

std::ostream& PointBounds::toStream(std::ostream& sl) const {
  detail::OstreamStateGuard guard{sl};
  sl << std::fixed << std::setprecision(7);
  sl << "Acts::PointBounds: (maxDistance) = ";
  sl << "(" << get(eR) << ")";
  return sl;
}

}  // namespace Acts
