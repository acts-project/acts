// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/EllipseBounds.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace Acts {

bool EllipseBounds::inside(const Vector2& lposition) const {
  double phi =
      detail::radian_sym(VectorHelpers::phi(lposition) - get(eAveragePhi));
  return (-get(eHalfPhiSector) <= phi) && (phi < get(eHalfPhiSector)) &&
         (square(lposition[eBoundLoc0] / get(eInnerRx)) +
          square(lposition[eBoundLoc1] / get(eInnerRy))) >= 1 &&
         (square(lposition[eBoundLoc0] / get(eOuterRx)) +
          square(lposition[eBoundLoc1] / get(eOuterRy))) < 1;
}

Vector2 EllipseBounds::closestPoint(
    const Vector2& /*lposition*/,
    const std::optional<SquareMatrix2>& /*metric*/) const {
  throw std::logic_error("Not implemented");
}

bool EllipseBounds::inside(const Vector2& lposition,
                           const BoundaryTolerance& boundaryTolerance) const {
  if (boundaryTolerance.isInfinite()) {
    return true;
  }

  if (auto absoluteBound = boundaryTolerance.asAbsoluteBoundOpt();
      absoluteBound.has_value()) {
    double tol0 = absoluteBound->tolerance0;
    double tol1 = absoluteBound->tolerance1;

    double phi =
        detail::radian_sym(VectorHelpers::phi(lposition) - get(eAveragePhi));
    double phiHalf = get(eHalfPhiSector) + tol1;

    bool insidePhi = (-phiHalf <= phi) && (phi < phiHalf);
    bool insideInner =
        (get(eInnerRx) <= tol0) || (get(eOuterRx) <= tol0) ||
        (1 < (square(lposition[eBoundLoc0] / (get(eInnerRx) - tol0)) +
              square(lposition[eBoundLoc1] / (get(eOuterRx) - tol0))));
    bool insideOuter =
        (square(lposition[eBoundLoc0] / (get(eInnerRy) + tol0)) +
         square(lposition[eBoundLoc1] / (get(eOuterRy) + tol0))) < 1;
    return insidePhi && insideInner && insideOuter;
  }

  throw std::logic_error("Unsupported boundary check type");
}

std::vector<Vector2> EllipseBounds::vertices(
    unsigned int quarterSegments) const {
  return detail::VerticesHelper::ellipsoidVertices(
      get(eInnerRx), get(eInnerRy), get(eOuterRx), get(eOuterRy),
      get(eAveragePhi), get(eHalfPhiSector), quarterSegments);
}

const RectangleBounds& EllipseBounds::boundingBox() const {
  return m_boundingBox;
}

std::ostream& EllipseBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::EllipseBounds:  (innerRadius0, outerRadius0, innerRadius1, "
        "outerRadius1, hPhiSector, averagePhi) = ";
  sl << "(" << get(eInnerRx) << ", " << get(eInnerRy) << ", " << get(eOuterRx)
     << ", " << get(eOuterRy) << ", " << get(eAveragePhi) << ", "
     << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ")";
  sl << std::setprecision(-1);
  return sl;
}

}  // namespace Acts
