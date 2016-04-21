///////////////////////////////////////////////////////////////////
// CylinderBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////

//Trk
#include "ACTS/Surfaces/CylinderBounds.h"
// STD/STL
#include <iostream>
#include <iomanip>
#include <math.h>

Acts::CylinderBounds::CylinderBounds() :
  m_boundValues(CylinderBounds::bv_length, 0.),
  m_checkPhi(false)
{}

Acts::CylinderBounds::CylinderBounds(double radius, double halez) :
    m_boundValues(CylinderBounds::bv_length, 0.),
    m_checkPhi(false)
{
    m_boundValues[CylinderBounds::bv_radius]        = fabs(radius);
    m_boundValues[CylinderBounds::bv_halfPhiSector] = M_PI;
    m_boundValues[CylinderBounds::bv_halfZ]         = fabs(halez);
}

Acts::CylinderBounds::CylinderBounds(double radius, double haphi, double halez) :
    m_boundValues(CylinderBounds::bv_length, 0.),
    m_checkPhi(true)
{
    m_boundValues[CylinderBounds::bv_radius]        = fabs(radius);
    m_boundValues[CylinderBounds::bv_halfPhiSector] = haphi;
    m_boundValues[CylinderBounds::bv_halfZ]         = fabs(halez);
}

Acts::CylinderBounds::CylinderBounds(double radius, double haphi, double averagephi, double halez) :
    m_boundValues(CylinderBounds::bv_length, 0.),
    m_checkPhi(true)
{
    m_boundValues[CylinderBounds::bv_radius]        = fabs(radius);
    m_boundValues[CylinderBounds::bv_averagePhi]    = averagephi;
    m_boundValues[CylinderBounds::bv_halfPhiSector] = haphi;
    m_boundValues[CylinderBounds::bv_halfZ]         = fabs(halez);
}

Acts::CylinderBounds::CylinderBounds(const Acts::CylinderBounds& cylbo) :
  Acts::SurfaceBounds(),
  m_boundValues(cylbo.m_boundValues),
  m_checkPhi(cylbo.m_checkPhi)
{}

Acts::CylinderBounds::~CylinderBounds()
{}

Acts::CylinderBounds& Acts::CylinderBounds::operator=(const Acts::CylinderBounds& cylbo)
{
  if (this!=&cylbo) {
    m_boundValues = cylbo.m_boundValues;
    m_checkPhi = cylbo.m_checkPhi;
  }
  return *this;
}

Acts::CylinderBounds& Acts::CylinderBounds::operator=(Acts::CylinderBounds&& cylbo)
{
  if (this!=&cylbo) {
    m_boundValues = std::move(cylbo.m_boundValues);
    m_checkPhi = cylbo.m_checkPhi;
  }
  return *this;
}

bool Acts::CylinderBounds::operator==(const SurfaceBounds& sbo) const
{
  // check the type first not to compare apples with oranges
  const Acts::CylinderBounds* cylbo = dynamic_cast<const Acts::CylinderBounds*>(&sbo);
  if (!cylbo) return false;
  return (m_boundValues == cylbo->m_boundValues);
}

double Acts::CylinderBounds::minDistance(const Acts::Vector2D& pos ) const
{
  const double pi2 = 2.*M_PI;

  double sZ = fabs(pos[Acts::eLOC_Z])-m_boundValues[CylinderBounds::bv_halfZ];
  double wF = m_boundValues[CylinderBounds::bv_halfPhiSector]; if(wF >= M_PI) return sZ;
  double dF = fabs(pos[Acts::eLOC_RPHI]/m_boundValues[CylinderBounds::bv_radius]-m_boundValues[CylinderBounds::bv_averagePhi]); if(dF>M_PI) dF=pi2-dF;
  double sF = 2.*m_boundValues[CylinderBounds::bv_radius]*sin(.5* (dF-wF));

  if(sF <= 0. || sZ <= 0.) {
      if (sF > sZ) return sF;
      else return sZ;
  }
  return sqrt(sF*sF+sZ*sZ);
}


// ostream operator overload
std::ostream& Acts::CylinderBounds::dump( std::ostream& sl ) const
{
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(7);
    sl << "Acts::CylinderBounds: (radius, averagePhi, halfPhiSector, halflengthInZ) = ";
    sl << "(" << this->r() << ", " << this->averagePhi() << ", ";
    sl << this->halfPhiSector() << ", " << this->halflengthZ() << ")";
    sl << std::setprecision(-1);
    return sl;
}






