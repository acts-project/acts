///////////////////////////////////////////////////////////////////
// DiscTrapezoidalBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////

//Trk
#include "Surfaces/DiscTrapezoidalBounds.h"
// STD/STL
#include <iostream>
#include <iomanip>

Acts::DiscTrapezoidalBounds::DiscTrapezoidalBounds() :
  Acts::DiscBounds(),
  m_boundValues(DiscTrapezoidalBounds::bv_length, 0.)
{}

Acts::DiscTrapezoidalBounds::DiscTrapezoidalBounds(double minhalfx, double maxhalfx, double maxR, double minR, double avephi, double stereo) :
  Acts::DiscBounds(),
  m_boundValues(DiscTrapezoidalBounds::bv_length, 0.)
{
  m_boundValues[DiscTrapezoidalBounds::bv_averagePhi]    = avephi;
  m_boundValues[DiscTrapezoidalBounds::bv_stereo]        = stereo;

  m_boundValues[DiscTrapezoidalBounds::bv_minHalfX]      = minhalfx;
  m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX]      = maxhalfx;
  if (m_boundValues[DiscTrapezoidalBounds::bv_minHalfX] >  m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX])
    swap(m_boundValues[DiscTrapezoidalBounds::bv_minHalfX],  m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX]);

  m_boundValues[DiscTrapezoidalBounds::bv_rMax] = minR;
  m_boundValues[DiscTrapezoidalBounds::bv_rMin] = maxR;
  if (m_boundValues[DiscTrapezoidalBounds::bv_rMin] >  m_boundValues[DiscTrapezoidalBounds::bv_rMax])
    swap(m_boundValues[DiscTrapezoidalBounds::bv_rMin],  m_boundValues[DiscTrapezoidalBounds::bv_rMax]);

  m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector] = asin(m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX]/m_boundValues[DiscTrapezoidalBounds::bv_rMax]);

  double hmax = sqrt(m_boundValues[DiscTrapezoidalBounds::bv_rMax]*m_boundValues[DiscTrapezoidalBounds::bv_rMax]-
		     m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX]*m_boundValues[DiscTrapezoidalBounds::bv_maxHalfX]);

  double hmin = sqrt(m_boundValues[DiscTrapezoidalBounds::bv_rMin]*m_boundValues[DiscTrapezoidalBounds::bv_rMin]-
		     m_boundValues[DiscTrapezoidalBounds::bv_minHalfX]*m_boundValues[DiscTrapezoidalBounds::bv_minHalfX]);

  m_boundValues[DiscTrapezoidalBounds::bv_rCenter] = (hmax+hmin)/2.;
  m_boundValues[DiscTrapezoidalBounds::bv_halfY]   = (hmax-hmin)/2.;


}

Acts::DiscTrapezoidalBounds::DiscTrapezoidalBounds(const DiscTrapezoidalBounds& disctrbo) :
  Acts::DiscBounds(),
  m_boundValues(disctrbo.m_boundValues)
{}

Acts::DiscTrapezoidalBounds::~DiscTrapezoidalBounds()
{}

Acts::DiscTrapezoidalBounds& Acts::DiscTrapezoidalBounds::operator=(const DiscTrapezoidalBounds& disctrbo)
{
  if (this!=&disctrbo)
    m_boundValues    = disctrbo.m_boundValues;
  return *this;
}

bool Acts::DiscTrapezoidalBounds::operator==(const Acts::SurfaceBounds& sbo) const
{
  // check the type first not to compare apples with oranges
  const Acts::DiscTrapezoidalBounds* disctrbo = dynamic_cast<const Acts::DiscTrapezoidalBounds*>(&sbo);
  if (!disctrbo) return false;
  return ( m_boundValues == disctrbo->m_boundValues);
}

double Acts::DiscTrapezoidalBounds::minDistance(const Acts::Vector2D& pos ) const
{
   const double pi2  = 2.*M_PI;
   double alpha = fabs(pos[Acts::eLOC_PHI]);
   if ( alpha>M_PI) alpha = pi2 - alpha;

   double r  = pos[Acts::eLOC_R];  if(r==0.) return m_boundValues[DiscTrapezoidalBounds::bv_rMin]*cos(m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector])/cos(alpha);

   // check if it is external in R
   double sr0 = m_boundValues[DiscTrapezoidalBounds::bv_rMin]*cos(m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector])/cos(alpha)-r;
   if(sr0 > 0.) return sr0;
   double sr1 = r-m_boundValues[DiscTrapezoidalBounds::bv_rMax]*cos(m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector])/cos(alpha);
   if(sr1 > 0.) return sr1;

   // check if it is external in phi
   if((alpha-m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector])>0.) return r*fabs(sin(alpha-m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector]));

   // if here it is inside:
   // Evaluate the distance from the 4 segments
   double dist [4] = { sr0,
		       sr1,
		       -r*sin(m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector]-alpha),
		       -r*sin(m_boundValues[DiscTrapezoidalBounds::bv_halfPhiSector]+alpha)};

   return *std::max_element(dist,dist+4);
}

// ostream operator overload
std::ostream& Acts::DiscTrapezoidalBounds::dump( std::ostream& sl ) const
{
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(7);
    sl << "Acts::DiscTrapezoidalBounds:  (innerRadius, outerRadius, hMinX, hMaxX, hlengthY, hPhiSector, averagePhi, rCenter, stereo) = ";
    sl << "(" << this->rMin() << ", " << this->rMax() << ", " << this->minHalflengthX() << ", " << this->maxHalflengthX() << ", " << this->halflengthY() << ", " << this->halfPhiSector() << ", " << this->averagePhi() << ", " << this->rCenter() << ", " << this->stereo() << ")";
    sl << std::setprecision(-1);
    return sl;
}
