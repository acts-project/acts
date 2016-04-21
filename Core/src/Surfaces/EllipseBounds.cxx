///////////////////////////////////////////////////////////////////
// EllipseBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/EllipseBounds.h"
// STD/STL
#include <iostream>
#include <iomanip>

Acts::EllipseBounds::EllipseBounds() :
    Acts::PlanarBounds(),
    m_boundValues(EllipseBounds::bv_length,0.)
{}


Acts::EllipseBounds::EllipseBounds(double minradX, double minradY, double maxradX, double maxradY, double hphisec) :
    Acts::PlanarBounds(),
    m_boundValues(EllipseBounds::bv_length,0.)
{
    m_boundValues[EllipseBounds::bv_rMinX]         = minradX;
    m_boundValues[EllipseBounds::bv_rMinY]         = minradY;
    m_boundValues[EllipseBounds::bv_rMaxX]         = maxradX;
    m_boundValues[EllipseBounds::bv_rMaxY]         = maxradY;
    m_boundValues[EllipseBounds::bv_averagePhi]    = 0.;
    m_boundValues[EllipseBounds::bv_halfPhiSector] = hphisec;
    if (m_boundValues[EllipseBounds::bv_rMinX] > m_boundValues[EllipseBounds::bv_rMaxX])
        swap(m_boundValues[EllipseBounds::bv_rMinX], m_boundValues[EllipseBounds::bv_rMaxX]);
    if (m_boundValues[EllipseBounds::bv_rMinY] > m_boundValues[EllipseBounds::bv_rMaxY])
        swap(m_boundValues[EllipseBounds::bv_rMinY], m_boundValues[EllipseBounds::bv_rMaxY]);
}

Acts::EllipseBounds::EllipseBounds(double minradX, double minradY, double maxradX, double maxradY, double avephi, double hphisec) :
    Acts::PlanarBounds(),
    m_boundValues(EllipseBounds::bv_length,0.)
{

    m_boundValues[EllipseBounds::bv_rMinX]         = minradX;
    m_boundValues[EllipseBounds::bv_rMinY]         = minradY;
    m_boundValues[EllipseBounds::bv_rMaxX]         = maxradX;
    m_boundValues[EllipseBounds::bv_rMaxY]         = maxradY;
    m_boundValues[EllipseBounds::bv_averagePhi]    = avephi;
    m_boundValues[EllipseBounds::bv_halfPhiSector] = hphisec;
    if (m_boundValues[EllipseBounds::bv_rMinX] > m_boundValues[EllipseBounds::bv_rMaxX])
        swap(m_boundValues[EllipseBounds::bv_rMinX], m_boundValues[EllipseBounds::bv_rMaxX]);
    if (m_boundValues[EllipseBounds::bv_rMinY] > m_boundValues[EllipseBounds::bv_rMaxY])
        swap(m_boundValues[EllipseBounds::bv_rMinY], m_boundValues[EllipseBounds::bv_rMaxY]);
}

Acts::EllipseBounds::EllipseBounds(const EllipseBounds& discbo) :
  Acts::PlanarBounds(),
  m_boundValues(discbo.m_boundValues)
{}

Acts::EllipseBounds::~EllipseBounds()
{}

Acts::EllipseBounds& Acts::EllipseBounds::operator=(const EllipseBounds& discbo)
{
  if (this!=&discbo)
      m_boundValues = discbo.m_boundValues;
  return *this;
}

Acts::EllipseBounds& Acts::EllipseBounds::operator=(EllipseBounds&& discbo)
{
  if (this!=&discbo)
    m_boundValues = std::move(discbo.m_boundValues);
  return *this;
}

bool Acts::EllipseBounds::operator==(const Acts::SurfaceBounds& sbo) const
{
  // check the type first not to compare apples with oranges
  const Acts::EllipseBounds* discbo = dynamic_cast<const Acts::EllipseBounds*>(&sbo);
  if (!discbo) return false;
  return (m_boundValues == discbo->m_boundValues);
}

// For ellipse bound this is only approximation which is valid
// only if m_boundValues[EllipseBounds::bv_rMinX] ~= m_boundValues[EllipseBounds::bv_rMinY]
// and m_boundValues[EllipseBounds::bv_rMaxX] ~= m_boundValues[EllipseBounds::bv_rMaxY]
//
double Acts::EllipseBounds::minDistance(const Acts::Vector2D& pos ) const
{
  const double pi2 = 2.*M_PI;

  double r    = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
  if(r==0.) {
      if (m_boundValues[EllipseBounds::bv_rMinX] <= m_boundValues[EllipseBounds::bv_rMinY]) return m_boundValues[EllipseBounds::bv_rMinX];
      return m_boundValues[EllipseBounds::bv_rMinY];
  }

  const double inv_r = 1. / r;
  double sn   = pos[1]*inv_r                     ;
  double cs   = pos[0]*inv_r                     ;
  double sf   = 0.                               ;
  double dF   = 0.                               ;

  if(m_boundValues[EllipseBounds::bv_halfPhiSector] < M_PI) {

    dF = atan2(cs,sn)-m_boundValues[EllipseBounds::bv_averagePhi];
    dF += (dF > M_PI) ? -pi2 : (dF < -M_PI) ? pi2 : 0;
    double df = fabs(dF)-m_boundValues[EllipseBounds::bv_halfPhiSector];
    sf = r*sin(df);
    if (df > 0.) r*=cos(df);
  }
  else {
    sf = -1.e+10;
  }

  if(sf <= 0. ) {

    double a   = cs/m_boundValues[EllipseBounds::bv_rMaxX]        ;
    double b   = sn/m_boundValues[EllipseBounds::bv_rMaxY]        ;
    double sr0 = r-1./sqrt(a*a+b*b); if(sr0 >=0.) return sr0;
    a          = cs/m_boundValues[EllipseBounds::bv_rMinX]        ;
    b          = sn/m_boundValues[EllipseBounds::bv_rMinY]        ;
    double sr1 = 1./sqrt(a*a+b*b)-r; if(sr1 >=0.) return sr1;
    if(sf < sr0) sf = sr0          ;
    if(sf < sr1) sf = sr1          ;
    return sf;
  }

  double fb;
  fb = (dF > 0.) ? m_boundValues[EllipseBounds::bv_averagePhi]+m_boundValues[EllipseBounds::bv_halfPhiSector] : m_boundValues[EllipseBounds::bv_averagePhi]-m_boundValues[EllipseBounds::bv_halfPhiSector];
  sn         = sin(fb)           ;
  cs         = cos(fb)           ;
  double a   = cs/m_boundValues[EllipseBounds::bv_rMaxX]        ;
  double b   = sn/m_boundValues[EllipseBounds::bv_rMaxY]        ;
  double sr0 = r-1./sqrt(a*a+b*b); if(sr0 >=0.) return sqrt(sr0*sr0+sf*sf);
  a          = cs/m_boundValues[EllipseBounds::bv_rMinX]        ;
  b          = sn/m_boundValues[EllipseBounds::bv_rMinY]        ;
  double sr1 = 1./sqrt(a*a+b*b)-r; if(sr1 >=0.) return sqrt(sr1*sr1+sf*sf);
  return sf;
}

// ostream operator overload
std::ostream& Acts::EllipseBounds::dump( std::ostream& sl ) const
{
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(7);
    sl << "Acts::EllipseBounds:  (innerRadiusX, innerRadiusY, outerRadiusX, outerRadiusY, hPhiSector) = ";
    sl << "(" << this->rMinX() << ", " << this->rMinY() << ", " << this->rMaxX() << ", " << this->rMaxY() << ", " << this->averagePhi() << ", "<< this->halfPhiSector() << ")";
    sl << std::setprecision(-1);
    return sl;
}

