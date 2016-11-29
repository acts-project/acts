// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RadialBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/RadialBounds.hpp"
#include <iomanip>
#include <iostream>
#include <cmath>

Acts::RadialBounds::RadialBounds(double minrad, double maxrad, double hphisec)
  : Acts::DiscBounds(RadialBounds::bv_length)
{
  m_valueStore.at(RadialBounds::bv_rMin)          = minrad;
  m_valueStore.at(RadialBounds::bv_rMax)          = maxrad;
  m_valueStore.at(RadialBounds::bv_averagePhi)    = 0.;
  m_valueStore.at(RadialBounds::bv_halfPhiSector) = hphisec;
  if (m_valueStore.at(RadialBounds::bv_rMin)
      > m_valueStore.at(RadialBounds::bv_rMax))
    std::swap(m_valueStore.at(RadialBounds::bv_rMin),
              m_valueStore.at(RadialBounds::bv_rMax));
}

Acts::RadialBounds::RadialBounds(double minrad,
                                 double maxrad,
                                 double avephi,
                                 double hphisec)
  : Acts::DiscBounds(RadialBounds::bv_length)
{
  m_valueStore.at(RadialBounds::bv_rMin)          = minrad;
  m_valueStore.at(RadialBounds::bv_rMax)          = maxrad;
  m_valueStore.at(RadialBounds::bv_averagePhi)    = avephi;
  m_valueStore.at(RadialBounds::bv_halfPhiSector) = hphisec;
  if (m_valueStore.at(RadialBounds::bv_rMin)
      > m_valueStore.at(RadialBounds::bv_rMax))
    std::swap(m_valueStore.at(RadialBounds::bv_rMin),
              m_valueStore.at(RadialBounds::bv_rMax));
}

Acts::RadialBounds::~RadialBounds()
{
}

Acts::RadialBounds&
Acts::RadialBounds::operator=(const RadialBounds& rbo)
{
  if (this != &rbo) SurfaceBounds::operator=(rbo);
  return *this;
}

double
Acts::RadialBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  const double pi2 = 2. * M_PI;

  double r = lpos[Acts::eLOC_R];
  if (r == 0.) return m_valueStore.at(RadialBounds::bv_rMin);
  double sf = 0.;
  double dF = 0.;

  if (m_valueStore.at(RadialBounds::bv_halfPhiSector) < M_PI) {
    dF = std::abs(lpos[Acts::eLOC_PHI]
              - m_valueStore.at(RadialBounds::bv_averagePhi));
    if (dF > M_PI) dF = pi2 - dF;
    dF -= m_valueStore.at(RadialBounds::bv_halfPhiSector);
    sf = r * sin(dF);
    if (dF > 0.) r *= cos(dF);

  } else {
    sf = -1.e+10;
  }

  if (sf <= 0.) {
    double sr0 = m_valueStore.at(RadialBounds::bv_rMin) - r;
    if (sr0 > 0.) return sr0;
    double sr1 = r - m_valueStore.at(RadialBounds::bv_rMax);
    if (sr1 > 0.) return sr1;
    if (sf < sr0) sf = sr0;
    if (sf < sr1) sf = sr1;
    return sf;
  }

  double sr0 = m_valueStore.at(RadialBounds::bv_rMin) - r;
  if (sr0 > 0.) return sqrt(sr0 * sr0 + sf * sf);
  double sr1 = r - m_valueStore.at(RadialBounds::bv_rMax);
  if (sr1 > 0.) return sqrt(sr1 * sr1 + sf * sf);
  return sf;
}

// ostream operator overload
std::ostream&
Acts::RadialBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::RadialBounds:  (innerRadius, outerRadius, hPhiSector) = ";
  sl << "(" << rMin() << ", " << rMax() << ", " << averagePhi() << ", "
     << halfPhiSector() << ")";
  sl << std::setprecision(-1);
  return sl;
}
