// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/DiamondBounds.hpp"

#include <iomanip>
#include <iostream>

Acts::SurfaceBounds::BoundsType Acts::DiamondBounds::type() const {
  return SurfaceBounds::eDiamond;
}

bool Acts::DiamondBounds::inside(const Acts::Vector2& lposition,
                                 const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(lposition, vertices());
}

std::vector<Acts::Vector2> Acts::DiamondBounds::vertices(
    unsigned int /*lseg*/) const {
  // Vertices starting at lower left (min rel. phi)
  // counter-clockwise
  double x1 = get(DiamondBounds::eHalfLengthXnegY);
  double y1 = get(DiamondBounds::eHalfLengthYneg);
  double x2 = get(DiamondBounds::eHalfLengthXzeroY);
  double y2 = 0.;
  double x3 = get(DiamondBounds::eHalfLengthXposY);
  double y3 = get(DiamondBounds::eHalfLengthYpos);
  return {{-x1, -y1}, {x1, -y1}, {x2, y2}, {x3, y3}, {-x3, y3}, {-x2, y2}};
}

const Acts::RectangleBounds& Acts::DiamondBounds::boundingBox() const {
  return m_boundingBox;
}

std::ostream& Acts::DiamondBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::DiamondBounds: (halfXatYneg, halfXatYzero, halfXatYpos, "
        "halfYneg, halfYpos) = ";
  sl << "(" << get(DiamondBounds::eHalfLengthXnegY) << ", "
     << get(DiamondBounds::eHalfLengthXzeroY) << ", "
     << get(DiamondBounds::eHalfLengthXposY) << ", "
     << get(DiamondBounds::eHalfLengthYneg) << ", "
     << get(DiamondBounds::eHalfLengthYpos) << ")";
  sl << std::setprecision(-1);
  return sl;
}
