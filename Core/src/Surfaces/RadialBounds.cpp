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
// STD/STL
#include <iostream>
#include <iomanip>

Acts::RadialBounds::RadialBounds() :
 Acts::DiscBounds(),
  m_boundValues(RadialBounds::bv_length, 0.)
{}

Acts::RadialBounds::RadialBounds(double minrad, double maxrad, double hphisec) :
  Acts::DiscBounds(),
  m_boundValues(RadialBounds::bv_length, 0.)
{
    m_boundValues.at(RadialBounds::bv_rMin)          = minrad;
    m_boundValues.at(RadialBounds::bv_rMax)          = maxrad;
    m_boundValues.at(RadialBounds::bv_averagePhi)    = 0.;
    m_boundValues.at(RadialBounds::bv_halfPhiSector) = hphisec;
    if (m_boundValues.at(RadialBounds::bv_rMin) >  m_boundValues.at(RadialBounds::bv_rMax))
        swap(m_boundValues.at(RadialBounds::bv_rMin),  m_boundValues.at(RadialBounds::bv_rMax));
}

Acts::RadialBounds::RadialBounds(double minrad, double maxrad, double avephi, double hphisec) :
  Acts::DiscBounds(),
  m_boundValues(RadialBounds::bv_length, 0.)
{
    m_boundValues.at(RadialBounds::bv_rMin)          = minrad;
    m_boundValues.at(RadialBounds::bv_rMax)          = maxrad;
    m_boundValues.at(RadialBounds::bv_averagePhi)    = avephi;
    m_boundValues.at(RadialBounds::bv_halfPhiSector) = hphisec;
    if (m_boundValues.at(RadialBounds::bv_rMin) >  m_boundValues.at(RadialBounds::bv_rMax))
        swap(m_boundValues.at(RadialBounds::bv_rMin),  m_boundValues.at(RadialBounds::bv_rMax));
}

Acts::RadialBounds::RadialBounds(const RadialBounds& discbo) :
  Acts::DiscBounds(),
  m_boundValues(discbo.m_boundValues)
{}

Acts::RadialBounds::~RadialBounds()
{}

Acts::RadialBounds& Acts::RadialBounds::operator=(const RadialBounds& discbo)
{
    if (this!=&discbo)
        m_boundValues    = discbo.m_boundValues;
    return *this;
}

bool Acts::RadialBounds::operator==(const Acts::SurfaceBounds& sbo) const
{
  // check the type first not to compare apples with oranges
  const Acts::RadialBounds* discbo = dynamic_cast<const Acts::RadialBounds*>(&sbo);
  if (!discbo) return false;
  return ( m_boundValues == discbo->m_boundValues);
}

double Acts::RadialBounds::minDistance(const Acts::Vector2D& pos ) const
{
  const double pi2  = 2.*M_PI;

  double r  = pos[Acts::eLOC_R];  if(r==0.) return m_boundValues.at(RadialBounds::bv_rMin);
  double sf = 0.       ;
  double dF = 0.       ;

  if(m_boundValues.at(RadialBounds::bv_halfPhiSector) < M_PI) {

    dF = fabs(pos[Acts::eLOC_PHI]-m_boundValues.at(RadialBounds::bv_averagePhi));
    if (dF>M_PI) dF=pi2-dF;
    dF-=m_boundValues.at(RadialBounds::bv_halfPhiSector);
    sf=r*sin(dF);
    if (dF > 0.) r*=cos(dF);

  }
  else {
    sf =-1.e+10;
  }

  if(sf <=0.) {

    double sr0 = m_boundValues.at(RadialBounds::bv_rMin)-r; if(sr0 > 0.) return sr0;
    double sr1 = r-m_boundValues.at(RadialBounds::bv_rMax); if(sr1 > 0.) return sr1;
    if(sf < sr0) sf = sr0;
    if(sf < sr1) sf = sr1;
    return sf;
  }

  double sr0 = m_boundValues.at(RadialBounds::bv_rMin)-r; if(sr0 > 0.) return sqrt(sr0*sr0+sf*sf);
  double sr1 = r-m_boundValues.at(RadialBounds::bv_rMax); if(sr1 > 0.) return sqrt(sr1*sr1+sf*sf);
  return sf;
}

// ostream operator overload
std::ostream& Acts::RadialBounds::dump( std::ostream& sl ) const
{
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(7);
    sl << "Acts::RadialBounds:  (innerRadius, outerRadius, hPhiSector) = ";
    sl << "(" << this->rMin() << ", " << this->rMax() << ", " << this->averagePhi() << ", "<< this->halfPhiSector() << ")";
    sl << std::setprecision(-1);
    return sl;
}

