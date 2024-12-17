// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RectangleBounds.hpp"

#include "Acts/Surfaces/detail/BoundaryCheckHelper.hpp"

#include <iomanip>
#include <iostream>

namespace Acts {

bool RectangleBounds::inside(const Vector2& lposition,
                             const BoundaryTolerance& boundaryTolerance) const {
  return detail::insideAlignedBox(m_min, m_max, boundaryTolerance, lposition,
                                  std::nullopt);
}

std::vector<Vector2> RectangleBounds::vertices(unsigned int /*lseg*/) const {
  // counter-clockwise starting from bottom-left corner
  return {m_min, {m_max.x(), m_min.y()}, m_max, {m_min.x(), m_max.y()}};
}

const RectangleBounds& RectangleBounds::boundingBox() const {
  return (*this);
}

// ostream operator overload
std::ostream& RectangleBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RectangleBounds:  (hlX, hlY) = "
     << "(" << 0.5 * (get(eMaxX) - get(eMinX)) << ", "
     << 0.5 * (get(eMaxY) - get(eMinY)) << ")";
  sl << "\n(lower left, upper right):\n";
  sl << min().transpose() << "\n" << max().transpose();
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
