// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::TrapezoidBounds::~TrapezoidBounds() = default;

Acts::SurfaceBounds::BoundsType Acts::TrapezoidBounds::type() const {
  return SurfaceBounds::eTrapezoid;
}

bool Acts::TrapezoidBounds::inside(const Acts::Vector2D& lposition,
                                   const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(lposition, vertices());
}

double Acts::TrapezoidBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(lposition, vertices());
}

std::vector<Acts::Vector2D> Acts::TrapezoidBounds::vertices(
    unsigned int /*lseg*/) const {
  double minhx = get(TrapezoidBounds::eHalfLengthXnegY);
  double maxhx = get(TrapezoidBounds::eHalfLengthXposY);
  double hy = get(TrapezoidBounds::eHalfLengthY);
  return {{-minhx, -hy}, {minhx, -hy}, {maxhx, hy}, {-maxhx, hy}};
}

const Acts::RectangleBounds& Acts::TrapezoidBounds::boundingBox() const {
  return m_boundingBox;
}

std::ostream& Acts::TrapezoidBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::TrapezoidBounds:  (halfXnegY, halfXposY, halfY) = "
     << "(" << get(eHalfLengthXnegY) << ", " << get(eHalfLengthXposY) << ", "
     << get(eHalfLengthY) << ")";
  sl << std::setprecision(-1);
  return sl;
}
