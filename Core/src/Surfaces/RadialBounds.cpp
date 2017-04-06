// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// RadialBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/RadialBounds.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "ACTS/Utilities/detail/periodic.hpp"

Acts::RadialBounds::RadialBounds(double minrad, double maxrad, double hphisec)
  : RadialBounds(minrad, maxrad, 0, hphisec)
{
}

Acts::RadialBounds::RadialBounds(double minrad,
                                 double maxrad,
                                 double avephi,
                                 double hphisec)
  : m_rMin(std::min(minrad, maxrad))
  , m_rMax(std::max(minrad, maxrad))
  , m_avgPhi(detail::radian_sym(avephi))
  , m_halfPhi(std::abs(hphisec))
{
}

Acts::RadialBounds::~RadialBounds()
{
}

Acts::RadialBounds*
Acts::RadialBounds::clone() const
{
  return new RadialBounds(*this);
}

Acts::SurfaceBounds::BoundsType
Acts::RadialBounds::type() const
{
  return SurfaceBounds::Disc;
}

std::vector<TDD_real_t>
Acts::RadialBounds::valueStore() const
{
  std::vector<TDD_real_t> values(RadialBounds::bv_length);
  values[RadialBounds::bv_rMin]          = rMin();
  values[RadialBounds::bv_rMax]          = rMax();
  values[RadialBounds::bv_averagePhi]    = averagePhi();
  values[RadialBounds::bv_halfPhiSector] = halfPhiSector();
  return values;
}

Acts::Vector2D
Acts::RadialBounds::shifted(const Acts::Vector2D& lpos) const
{
  Vector2D tmp;
  tmp[eLOC_R]   = lpos[eLOC_R];
  tmp[eLOC_PHI] = detail::radian_sym(lpos[eLOC_PHI] - averagePhi());
  return tmp;
}

bool
Acts::RadialBounds::inside(const Acts::Vector2D&      lpos,
                           const Acts::BoundaryCheck& bcheck) const
{
  return bcheck.isInside(
      shifted(lpos), rMin(), rMax(), -halfPhiSector(), halfPhiSector());
}

double
Acts::RadialBounds::distanceToBoundary(const Acts::Vector2D& lpos) const
{
  return BoundaryCheck(true).distance(
      shifted(lpos), rMin(), rMax(), -halfPhiSector(), halfPhiSector());
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
