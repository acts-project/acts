// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::SurfaceBounds::BoundsType Acts::RadialBounds::type() const {
  return SurfaceBounds::eDisc;
}

Acts::Vector2D Acts::RadialBounds::shifted(
    const Acts::Vector2D& lposition) const {
  Vector2D tmp;
  tmp[eLOC_R] = lposition[eLOC_R];
  tmp[eLOC_PHI] = detail::radian_sym(lposition[eLOC_PHI] - get(eAveragePhi));
  return tmp;
}

bool Acts::RadialBounds::inside(const Acts::Vector2D& lposition,
                                const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(shifted(lposition),
                         Vector2D(get(eMinR), -get(eHalfPhiSector)),
                         Vector2D(get(eMaxR), get(eHalfPhiSector)));
}

double Acts::RadialBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(
      shifted(lposition), Vector2D(get(eMinR), -get(eHalfPhiSector)),
      Vector2D(get(eMaxR), get(eHalfPhiSector)));
}

std::vector<Acts::Vector2D> Acts::RadialBounds::vertices(
    unsigned int lseg) const {
  return detail::VerticesHelper::circularVertices(
      get(eMinR), get(eMaxR), get(eAveragePhi), get(eHalfPhiSector), lseg);
}

std::ostream& Acts::RadialBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RadialBounds:  (innerRadius, outerRadius, hPhiSector, "
        "averagePhi) = ";
  sl << "(" << get(eMinR) << ", " << get(eMaxR) << ", " << get(eHalfPhiSector)
     << ", " << get(eAveragePhi) << ")";
  sl << std::setprecision(-1);
  return sl;
}
