// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/CylinderBounds.hpp"
#include <iomanip>
#include <iostream>
#include <math.h>

Acts::CylinderBounds::CylinderBounds(double radius, double halez)
  : m_valueStore(CylinderBounds::bv_length, 0.), m_checkPhi(false)
{
  m_valueStore.at(CylinderBounds::bv_radius)        = fabs(radius);
  m_valueStore.at(CylinderBounds::bv_halfPhiSector) = M_PI;
  m_valueStore.at(CylinderBounds::bv_halfZ)         = fabs(halez);
}

Acts::CylinderBounds::CylinderBounds(double radius, double haphi, double halez)
  : m_valueStore(CylinderBounds::bv_length, 0.), m_checkPhi(true)
{
  m_valueStore.at(CylinderBounds::bv_radius)        = fabs(radius);
  m_valueStore.at(CylinderBounds::bv_halfPhiSector) = haphi;
  m_valueStore.at(CylinderBounds::bv_halfZ)         = fabs(halez);
}

Acts::CylinderBounds::CylinderBounds(double radius,
                                     double haphi,
                                     double averagephi,
                                     double halez)
  : m_valueStore(CylinderBounds::bv_length, 0.), m_checkPhi(true)
{
  m_valueStore.at(CylinderBounds::bv_radius)        = fabs(radius);
  m_valueStore.at(CylinderBounds::bv_averagePhi)    = averagephi;
  m_valueStore.at(CylinderBounds::bv_halfPhiSector) = haphi;
  m_valueStore.at(CylinderBounds::bv_halfZ)         = fabs(halez);
}

Acts::CylinderBounds::CylinderBounds(const Acts::CylinderBounds& cylbo)
  : Acts::SurfaceBounds()
  , m_valueStore(cylbo.m_valueStore)
  , m_checkPhi(cylbo.m_checkPhi)
{
}

Acts::CylinderBounds::~CylinderBounds()
{
}

Acts::CylinderBounds&
Acts::CylinderBounds::operator=(const Acts::CylinderBounds& cylbo)
{
  if (this != &cylbo) {
    m_valueStore = cylbo.m_valueStore;
    m_checkPhi    = cylbo.m_checkPhi;
  }
  return *this;
}

double
Acts::CylinderBounds::minDistance(const Acts::Vector2D& pos) const
{
  const double pi2 = 2. * M_PI;

  double sZ
      = fabs(pos[Acts::eLOC_Z]) - m_valueStore.at(CylinderBounds::bv_halfZ);
  double wF = m_valueStore.at(CylinderBounds::bv_halfPhiSector);
  if (wF >= M_PI) return sZ;
  double dF
      = fabs(pos[Acts::eLOC_RPHI] / m_valueStore.at(CylinderBounds::bv_radius)
             - m_valueStore.at(CylinderBounds::bv_averagePhi));
  if (dF > M_PI) dF = pi2 - dF;
  double sF
      = 2. * m_valueStore.at(CylinderBounds::bv_radius) * sin(.5 * (dF - wF));

  if (sF <= 0. || sZ <= 0.) {
    if (sF > sZ)
      return sF;
    else
      return sZ;
  }
  return sqrt(sF * sF + sZ * sZ);
}

// ostream operator overload
std::ostream&
Acts::CylinderBounds::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(7);
  sl << "Acts::CylinderBounds: (radius, averagePhi, halfPhiSector, "
        "halflengthInZ) = ";
  sl << "(" << this->r() << ", " << this->averagePhi() << ", ";
  sl << this->halfPhiSector() << ", " << this->halflengthZ() << ")";
  sl << std::setprecision(-1);
  return sl;
}
