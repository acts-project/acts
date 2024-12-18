// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/RadialBounds.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/detail/BoundaryCheckHelper.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <iomanip>
#include <iostream>

namespace Acts {

SquareMatrix2 RadialBounds::boundToCartesianJacobian(
    const Vector2& /*lposition*/) const {
  return SquareMatrix2::Identity();  // TODO
}

SquareMatrix2 RadialBounds::cartesianToBoundJacobian(
    const Vector2& /*lposition*/) const {
  return SquareMatrix2::Identity();  // TODO
}

Vector2 RadialBounds::shifted(const Vector2& lposition) const {
  Vector2 tmp;
  tmp[eBoundLoc0] = lposition[eBoundLoc0];
  tmp[eBoundLoc1] =
      detail::radian_sym(lposition[eBoundLoc1] - get(eAveragePhi));
  return tmp;
}

bool RadialBounds::inside(const Vector2& lposition) const {
  return detail::insideAlignedBox(Vector2(get(eMinR), -get(eHalfPhiSector)),
                                  Vector2(get(eMaxR), get(eHalfPhiSector)),
                                  BoundaryTolerance::None(), shifted(lposition),
                                  std::nullopt);
}

Vector2 RadialBounds::closestPoint(
    const Vector2& lposition,
    const std::optional<SquareMatrix2>& metric) const {
  return detail::computeClosestPointOnAlignedBox(
      Vector2(get(eMinR), -get(eHalfPhiSector)),
      Vector2(get(eMaxR), get(eHalfPhiSector)), shifted(lposition), metric);
}

bool RadialBounds::inside(const Vector2& lposition,
                          const BoundaryTolerance& boundaryTolerance) const {
  return detail::insideAlignedBox(Vector2(get(eMinR), -get(eHalfPhiSector)),
                                  Vector2(get(eMaxR), get(eHalfPhiSector)),
                                  boundaryTolerance, shifted(lposition),
                                  std::nullopt);
}

std::vector<Vector2> RadialBounds::vertices(unsigned int lseg) const {
  return detail::VerticesHelper::circularVertices(
      get(eMinR), get(eMaxR), get(eAveragePhi), get(eHalfPhiSector), lseg);
}

std::ostream& RadialBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "RadialBounds:  (innerRadius, outerRadius, hPhiSector, "
        "averagePhi) = ";
  sl << "(" << get(eMinR) << ", " << get(eMaxR) << ", " << get(eHalfPhiSector)
     << ", " << get(eAveragePhi) << ")";
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
