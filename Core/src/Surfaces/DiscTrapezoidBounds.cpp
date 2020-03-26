// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/DiscTrapezoidBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

Acts::DiscTrapezoidBounds::DiscTrapezoidBounds(double halfXminR,
                                               double halfXmaxR, double minR,
                                               double maxR, double avgPhi,
                                               double stereo) noexcept(false)
    : m_values({halfXminR, halfXmaxR, minR, maxR, avgPhi, stereo}) {
  checkConsistency();
}

Acts::SurfaceBounds::BoundsType Acts::DiscTrapezoidBounds::type() const {
  return SurfaceBounds::eDiscTrapezoid;
}

Acts::Vector2D Acts::DiscTrapezoidBounds::toLocalCartesian(
    const Acts::Vector2D& lposition) const {
  return {lposition[eLOC_R] * std::sin(lposition[eLOC_PHI] - get(eAveragePhi)),
          lposition[eLOC_R] * std::cos(lposition[eLOC_PHI] - get(eAveragePhi))};
}

Acts::ActsMatrixD<2, 2> Acts::DiscTrapezoidBounds::jacobianToLocalCartesian(
    const Acts::Vector2D& lposition) const {
  ActsMatrixD<2, 2> jacobian;
  jacobian(0, eLOC_R) = std::sin(lposition[eLOC_PHI] - get(eAveragePhi));
  jacobian(1, eLOC_R) = std::cos(lposition[eLOC_PHI] - get(eAveragePhi));
  jacobian(0, eLOC_PHI) = lposition[eLOC_R] * std::cos(lposition[eLOC_PHI]);
  jacobian(1, eLOC_PHI) = lposition[eLOC_R] * -std::sin(lposition[eLOC_PHI]);
  return jacobian;
}

bool Acts::DiscTrapezoidBounds::inside(
    const Acts::Vector2D& lposition, const Acts::BoundaryCheck& bcheck) const {
  Vector2D vertices[] = {{get(eHalfLengthXminR), get(eMinR)},
                         {get(eHalfLengthXmaxR), get(eMaxR)},
                         {-get(eHalfLengthXmaxR), get(eMaxR)},
                         {-get(eHalfLengthXminR), get(eMinR)}};
  auto jacobian = jacobianToLocalCartesian(lposition);
  return bcheck.transformed(jacobian).isInside(toLocalCartesian(lposition),
                                               vertices);
}

double Acts::DiscTrapezoidBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  Vector2D vertices[] = {{get(eHalfLengthXminR), get(eMinR)},
                         {get(eHalfLengthXmaxR), get(eMaxR)},
                         {-get(eHalfLengthXmaxR), get(eMaxR)},
                         {-get(eHalfLengthXminR), get(eMinR)}};
  return BoundaryCheck(true).distance(toLocalCartesian(lposition), vertices);
}

std::vector<Acts::Vector2D> Acts::DiscTrapezoidBounds::vertices(
    unsigned int /*lseg*/) const {
  Vector2D cAxis(std::cos(get(eAveragePhi)), std::sin(get(eAveragePhi)));
  Vector2D nAxis(cAxis.y(), -cAxis.x());
  return {get(eMinR) * cAxis - get(eHalfLengthXminR) * nAxis,
          get(eMinR) * cAxis + get(eHalfLengthXminR) * nAxis,
          get(eMaxR) * cAxis + get(eHalfLengthXmaxR) * nAxis,
          get(eMaxR) * cAxis - get(eHalfLengthXmaxR) * nAxis};
}

// ostream operator overload
std::ostream& Acts::DiscTrapezoidBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::DiscTrapezoidBounds: (innerRadius, outerRadius, "
        "halfLengthXminR, "
        "halfLengthXmaxR, halfLengthY, halfPhiSector, averagePhi, rCenter, "
        "stereo) = ";
  sl << "(" << get(eMinR) << ", " << get(eMaxR) << ", " << get(eHalfLengthXminR)
     << ", " << get(eHalfLengthXmaxR) << ", " << halfLengthY() << ", "
     << halfPhiSector() << ", " << get(eAveragePhi) << ", " << rCenter() << ", "
     << stereo() << ")";
  sl << std::setprecision(-1);
  return sl;
}
