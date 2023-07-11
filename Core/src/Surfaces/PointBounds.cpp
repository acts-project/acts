// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/PointBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::SurfaceBounds::BoundsType Acts::PointBounds::type() const {
  return SurfaceBounds::eLine;
}

bool Acts::PointBounds::inside(const Acts::Vector2& lposition,
                               const Acts::BoundaryCheck& bcheck) const {
  double r = get(PointBounds::eR);
  return bcheck.isInside(lposition, Vector2(0, 0), Vector2(r, 2 * M_PI));
}

// ostream operator overload
std::ostream& Acts::PointBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::PointBounds: (radius) = ";
  sl << "(" << get(PointBounds::eR) << ")";
  sl << std::setprecision(-1);
  return sl;
}
