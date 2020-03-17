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

Acts::RadialBounds::RadialBounds(double minrad, double maxrad, double hphisec)
    : RadialBounds(minrad, maxrad, 0, hphisec) {}

Acts::RadialBounds::RadialBounds(double minrad, double maxrad, double avephi,
                                 double hphisec)
    : m_rMin(std::min(minrad, maxrad)),
      m_rMax(std::max(minrad, maxrad)),
      m_avgPhi(detail::radian_sym(avephi)),
      m_halfPhi(std::abs(hphisec)) {}

Acts::RadialBounds* Acts::RadialBounds::clone() const {
  return new RadialBounds(*this);
}

Acts::SurfaceBounds::BoundsType Acts::RadialBounds::type() const {
  return SurfaceBounds::Disc;
}

std::vector<TDD_real_t> Acts::RadialBounds::valueStore() const {
  std::vector<TDD_real_t> values(RadialBounds::bv_length);
  values[RadialBounds::bv_rMin] = rMin();
  values[RadialBounds::bv_rMax] = rMax();
  values[RadialBounds::bv_averagePhi] = averagePhi();
  values[RadialBounds::bv_halfPhiSector] = halfPhiSector();
  return values;
}

Acts::Vector2D Acts::RadialBounds::shifted(
    const Acts::Vector2D& lposition) const {
  Vector2D tmp;
  tmp[eLOC_R] = lposition[eLOC_R];
  tmp[eLOC_PHI] = detail::radian_sym(lposition[eLOC_PHI] - averagePhi());
  return tmp;
}

bool Acts::RadialBounds::inside(const Acts::Vector2D& lposition,
                                const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(shifted(lposition), Vector2D(rMin(), -halfPhiSector()),
                         Vector2D(rMax(), halfPhiSector()));
}

double Acts::RadialBounds::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(shifted(lposition),
                                      Vector2D(rMin(), -halfPhiSector()),
                                      Vector2D(rMax(), halfPhiSector()));
}

std::vector<Acts::Vector2D> Acts::RadialBounds::vertices(
    unsigned int lseg) const {
  return detail::VerticesHelper::circularVertices(m_rMin, m_rMax, m_avgPhi,
                                                  m_halfPhi, lseg);
}

// ostream operator overload
std::ostream& Acts::RadialBounds::toStream(std::ostream& sl) const {
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RadialBounds:  (innerRadius, outerRadius, hPhiSector) = ";
  sl << "(" << rMin() << ", " << rMax() << ", " << averagePhi() << ", "
     << halfPhiSector() << ")";
  sl << std::setprecision(-1);
  return sl;
}
