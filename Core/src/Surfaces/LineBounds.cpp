// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/LineBounds.hpp"

#include <iomanip>
#include <iostream>

Acts::SurfaceBounds::BoundsType Acts::LineBounds::type() const {
  return SurfaceBounds::eLine;
}

bool Acts::LineBounds::inside(const Acts::Vector2D& lposition,
                              const Acts::BoundaryCheck& bcheck) const {
  double r = get(LineBounds::eR);
  double halfLengthZ = get(LineBounds::eHalfLengthZ);
  return bcheck.isInside(lposition, Vector2D(-r, -halfLengthZ),
                         Vector2D(r, halfLengthZ));
}

double Acts::LineBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  // per definition the min Distance of a correct local position is r
  return lposition[Acts::eLOC_R];
}

// ostream operator overload
std::ostream& Acts::LineBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::LineBounds: (radius, halflengthInZ) = ";
  sl << "(" << get(LineBounds::eR) << ", " << get(LineBounds::eHalfLengthZ)
     << ")";
  sl << std::setprecision(-1);
  return sl;
}
