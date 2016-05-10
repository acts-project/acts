///////////////////////////////////////////////////////////////////
// DiamondBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/DiamondBounds.h"
// STD/STL
#include <iostream>
#include <iomanip>
#include <math.h>

// default constructor
Acts::DiamondBounds::DiamondBounds() :
    Acts::PlanarBounds(),
    m_boundValues(DiamondBounds::bv_length, 0.),
    m_alpha1(0.),
    m_alpha2(0.)
{}

// constructor from arguments I
Acts::DiamondBounds::DiamondBounds(double minhalex, double medhalex, double maxhalex, double haley1, double haley2) :
    Acts::PlanarBounds(),
    m_boundValues(DiamondBounds::bv_length, 0.),
    m_alpha1(0.),
    m_alpha2(0.)
{
  m_boundValues[DiamondBounds::bv_minHalfX] = minhalex;
  m_boundValues[DiamondBounds::bv_medHalfX] = medhalex;
  m_boundValues[DiamondBounds::bv_maxHalfX] = maxhalex;
  m_boundValues[DiamondBounds::bv_halfY1]   = haley1;
  m_boundValues[DiamondBounds::bv_halfY2]   = haley2;
  if (minhalex > maxhalex)
      swap(m_boundValues[DiamondBounds::bv_minHalfX], m_boundValues[DiamondBounds::bv_maxHalfX]);
  initCache();

}

// copy constructor
Acts::DiamondBounds::DiamondBounds(const DiamondBounds& diabo) :
  Acts::PlanarBounds(),
  m_boundValues(diabo.m_boundValues),
  m_alpha1(diabo.m_alpha1),
  m_alpha2(diabo.m_alpha2)
{}

// destructor
Acts::DiamondBounds::~DiamondBounds()
{}

Acts::DiamondBounds& Acts::DiamondBounds::operator=(const DiamondBounds& diabo)
{
  if (this!=&diabo){
      m_boundValues = diabo.m_boundValues;
      m_alpha1      = diabo.m_alpha1;
      m_alpha2      = diabo.m_alpha2;
  }
  return *this;
}

bool Acts::DiamondBounds::operator==(const Acts::SurfaceBounds& sbo) const
{
  // check the type first not to compare apples with oranges
  const Acts::DiamondBounds* diabo = dynamic_cast<const Acts::DiamondBounds*>(&sbo);
  if (!diabo) return false;
  return (m_boundValues == diabo->m_boundValues);
  }

// checking if inside bounds (Full symmetrical Diamond)
bool Acts::DiamondBounds::inside(const Acts::Vector2D& locpo, double tol1, double tol2) const
{
    return this->insideFull(locpo, tol1, tol2);
}

// checking if inside bounds (Full symmetrical Diamond)
bool Acts::DiamondBounds::insideFull(const Acts::Vector2D& locpo, double tol1, double tol2) const
{
  // the cases:
  // (0)
  if (!m_boundValues[DiamondBounds::bv_halfY1] && !m_boundValues[DiamondBounds::bv_minHalfX]) return false;
  // (1)
  if (locpo[Acts::eLOC_Y] < -2. * m_boundValues[DiamondBounds::bv_halfY1] - tol2) return false;
  if (locpo[Acts::eLOC_Y] >  2. * m_boundValues[DiamondBounds::bv_halfY2] + tol2) return false;
  // (2)
  if (fabs(locpo[Acts::eLOC_X]) > (m_boundValues[DiamondBounds::bv_medHalfX] + tol1)) return false;
  // (3)
  if (fabs(locpo[Acts::eLOC_X]) < (fmin(m_boundValues[DiamondBounds::bv_minHalfX],m_boundValues[DiamondBounds::bv_maxHalfX]) - tol1)) return true;
  // (4)
  /** use basic calculation of a straight line */
  if (locpo[Acts::eLOC_Y]<0) {
    double k = m_boundValues[DiamondBounds::bv_halfY1]>0. ? (m_boundValues[DiamondBounds::bv_medHalfX] - m_boundValues[DiamondBounds::bv_minHalfX])/2/m_boundValues[DiamondBounds::bv_halfY1] : 0.;
    return ( fabs(locpo[Acts::eLOC_X]) <= m_boundValues[DiamondBounds::bv_medHalfX]-k*fabs(locpo[Acts::eLOC_Y]) );
  } else {
    double k = m_boundValues[DiamondBounds::bv_halfY2]>0. ? (m_boundValues[DiamondBounds::bv_medHalfX] - m_boundValues[DiamondBounds::bv_maxHalfX])/2/m_boundValues[DiamondBounds::bv_halfY2] : 0.;
    return ( fabs(locpo[Acts::eLOC_X]) <= m_boundValues[DiamondBounds::bv_medHalfX]-k*fabs(locpo[Acts::eLOC_Y]) );
  }
}

// opening angle in point A
double Acts::DiamondBounds::alpha1() const
{
   return m_alpha1;
}

// opening angle in point A'
double Acts::DiamondBounds::alpha2() const
{
   return m_alpha2;
}

double Acts::DiamondBounds::minDistance(const Acts::Vector2D& pos ) const
{
  const int Np = 6;

  double y1 = 2.*m_boundValues[DiamondBounds::bv_halfY1];
  double y2 = 2.*m_boundValues[DiamondBounds::bv_halfY2];

  double X [6] = {-m_boundValues[DiamondBounds::bv_minHalfX],
                  -m_boundValues[DiamondBounds::bv_medHalfX],
                  -m_boundValues[DiamondBounds::bv_maxHalfX],
                  m_boundValues[DiamondBounds::bv_maxHalfX],
                  m_boundValues[DiamondBounds::bv_medHalfX],
                  m_boundValues[DiamondBounds::bv_minHalfX]};
  double Y [6] = {        -y1,         0.,         y2,        y2,        0.,       -y1};

  double dm = 1.e+20;
  double Ao =     0.;
  bool   in =   true;

  for (int i=0; i!=Np; ++i) {

    int j = i+1; if (j==Np) j=0;

    double x  = X[i]-pos[0];
    double y  = Y[i]-pos[1];
    double dx = X[j]-X[i]  ;
    double dy = Y[j]-Y[i]  ;
    double A  =  x*dy-y*dx ;
    double S  =-(x*dx+y*dy);

    if (S <= 0.) {
        double d = x*x+y*y;
        if (d<dm) dm=d;
    }
    else {
      double a = dx*dx+dy*dy;
      if (S <= a) {
          double d = (A*A)/a;
          if (d<dm) dm=d;
      }
    }
    if (i && in && Ao*A < 0.) in = false;
    Ao = A;
  }
  if (in) return -sqrt(dm);
  return sqrt(dm);
}

// ostream operator overload
std::ostream& Acts::DiamondBounds::dump( std::ostream& sl ) const
{
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(7);
    sl << "Acts::DiamondBounds:  (minHlenghtX, medHlengthX, maxHlengthX, hlengthY1, hlengthY2 ) = ";
    sl << "(" << m_boundValues[DiamondBounds::bv_minHalfX] << ", "
        << m_boundValues[DiamondBounds::bv_medHalfX] <<", "
        << m_boundValues[DiamondBounds::bv_maxHalfX] << ", "
        << m_boundValues[DiamondBounds::bv_halfY1] << ", "
        << m_boundValues[DiamondBounds::bv_halfY2] << ")";
    sl << std::setprecision(-1);
    return sl;
}

