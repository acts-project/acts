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

#include "Acts/Utilities/detail/periodic.hpp"

Acts::DiscTrapezoidBounds::DiscTrapezoidBounds(double minhalfx, double maxhalfx,
                                               double minR, double maxR,
                                               double avephi, double stereo)
    : m_rMin(std::min(std::abs(minR), std::abs(maxR))),
      m_rMax(std::max(std::abs(minR), std::abs(maxR))),
      m_minHalfX(std::abs(minhalfx)),
      m_maxHalfX(std::abs(maxhalfx)),
      m_avgPhi(detail::radian_sym(avephi)),
      m_stereo(stereo) {}

Acts::DiscTrapezoidBounds* Acts::DiscTrapezoidBounds::clone() const {
  return new DiscTrapezoidBounds(*this);
}

Acts::SurfaceBounds::BoundsType Acts::DiscTrapezoidBounds::type() const {
  return SurfaceBounds::DiscTrapezoidal;
}

std::vector<TDD_real_t> Acts::DiscTrapezoidBounds::valueStore() const {
  std::vector<TDD_real_t> values(DiscTrapezoidBounds::bv_length);
  values[bv_rMin] = rMin();
  values[bv_rMax] = rMax();
  values[bv_minHalfX] = minHalflengthX();
  values[bv_maxHalfX] = maxHalflengthX();
  values[bv_averagePhi] = averagePhi();
  values[bv_stereo] = m_stereo;
  return values;
}

Acts::Vector2D Acts::DiscTrapezoidBounds::toLocalCartesian(
    const Acts::Vector2D& lposition) const {
  return {lposition[eLOC_R] * std::sin(lposition[eLOC_PHI] - m_avgPhi),
          lposition[eLOC_R] * std::cos(lposition[eLOC_PHI] - m_avgPhi)};
}

Acts::ActsMatrixD<2, 2> Acts::DiscTrapezoidBounds::jacobianToLocalCartesian(
    const Acts::Vector2D& lposition) const {
  ActsMatrixD<2, 2> jacobian;
  jacobian(0, eLOC_R) = std::sin(lposition[eLOC_PHI] - m_avgPhi);
  jacobian(1, eLOC_R) = std::cos(lposition[eLOC_PHI] - m_avgPhi);
  jacobian(0, eLOC_PHI) = lposition[eLOC_R] * std::cos(lposition[eLOC_PHI]);
  jacobian(1, eLOC_PHI) = lposition[eLOC_R] * -std::sin(lposition[eLOC_PHI]);
  return jacobian;
}

bool Acts::DiscTrapezoidBounds::inside(
    const Acts::Vector2D& lposition, const Acts::BoundaryCheck& bcheck) const {
  Vector2D vertices[] = {{minHalflengthX(), rMin()},
                         {maxHalflengthX(), rMax()},
                         {-maxHalflengthX(), rMax()},
                         {-minHalflengthX(), rMin()}};
  auto jacobian = jacobianToLocalCartesian(lposition);
  return bcheck.transformed(jacobian).isInside(toLocalCartesian(lposition),
                                               vertices);
}

double Acts::DiscTrapezoidBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  Vector2D vertices[] = {{minHalflengthX(), rMin()},
                         {maxHalflengthX(), rMax()},
                         {-maxHalflengthX(), rMax()},
                         {-minHalflengthX(), rMin()}};
  return BoundaryCheck(true).distance(toLocalCartesian(lposition), vertices);
}

std::vector<Acts::Vector2D> Acts::DiscTrapezoidBounds::vertices(
    unsigned int /*lseg*/) const {
  Vector2D cAxis(std::cos(m_avgPhi), std::sin(m_avgPhi));
  Vector2D nAxis(cAxis.y(), -cAxis.x());
  return {
      m_rMin * cAxis - m_minHalfX * nAxis, m_rMin * cAxis + m_minHalfX * nAxis,
      m_rMax * cAxis + m_maxHalfX * nAxis, m_rMax * cAxis - m_maxHalfX * nAxis};
}

// ostream operator overload
std::ostream& Acts::DiscTrapezoidBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::DiscTrapezoidBounds:  (innerRadius, outerRadius, hMinX, "
        "hMaxX, hlengthY, hPhiSector, averagePhi, rCenter, stereo) = ";
  sl << "(" << rMin() << ", " << rMax() << ", " << minHalflengthX() << ", "
     << maxHalflengthX() << ", " << halflengthY() << ", " << halfPhiSector()
     << ", " << averagePhi() << ", " << rCenter() << ", " << stereo() << ")";
  sl << std::setprecision(-1);
  return sl;
}
