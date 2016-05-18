///////////////////////////////////////////////////////////////////
// TrapezoidBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/TrapezoidBounds.hpp"
// STD/STL
#include <iostream>
#include <iomanip>
#include <math.h>

// default constructor
Acts::TrapezoidBounds::TrapezoidBounds() :
    Acts::PlanarBounds(),
    m_boundValues(TrapezoidBounds::bv_length, 0.),
    m_alpha(0.),
    m_beta(0.)
{}

// constructor from arguments I
Acts::TrapezoidBounds::TrapezoidBounds(double minhalex, double maxhalex, double haley) :
    Acts::PlanarBounds(),
    m_boundValues(TrapezoidBounds::bv_length, 0.),
    m_alpha(0.),
    m_beta(0.)
{
    m_boundValues.at(TrapezoidBounds::bv_minHalfX) = fabs(minhalex);
    m_boundValues.at(TrapezoidBounds::bv_maxHalfX) = fabs(maxhalex);
    m_boundValues.at(TrapezoidBounds::bv_halfY)   = fabs(haley);
    if (m_boundValues.at(TrapezoidBounds::bv_minHalfX) > m_boundValues.at(TrapezoidBounds::bv_maxHalfX))
        swap(m_boundValues.at(TrapezoidBounds::bv_minHalfX), m_boundValues.at(TrapezoidBounds::bv_maxHalfX));

}

// constructor from arguments II
Acts::TrapezoidBounds::TrapezoidBounds(double minhalex, double haley, double alpha, double beta) :
    Acts::PlanarBounds(),
    m_boundValues(TrapezoidBounds::bv_length, 0.),
    m_alpha(alpha),
    m_beta(beta)
{
    double gamma = (alpha > beta) ? (alpha - 0.5*M_PI) : (beta - 0.5*M_PI);
    // now fill them
    m_boundValues.at(TrapezoidBounds::bv_minHalfX) = fabs(minhalex);
    m_boundValues.at(TrapezoidBounds::bv_maxHalfX) = minhalex + (2.*m_boundValues.at(TrapezoidBounds::bv_halfY))*tan(gamma);
    m_boundValues.at(TrapezoidBounds::bv_halfY)   = fabs(haley);
}

// copy constructor
Acts::TrapezoidBounds::TrapezoidBounds(const TrapezoidBounds& trabo) :
    Acts::PlanarBounds(),
    m_boundValues(trabo.m_boundValues),
    m_alpha(trabo.m_alpha),
    m_beta(trabo.m_beta)
{}

// destructor
Acts::TrapezoidBounds::~TrapezoidBounds()
{}

Acts::TrapezoidBounds& Acts::TrapezoidBounds::operator=(const TrapezoidBounds& trabo)
{
    if (this!=&trabo){
        m_boundValues = trabo.m_boundValues;
        m_alpha       = trabo.m_alpha;
        m_beta        = trabo.m_beta;
    }
    return *this;
}

bool Acts::TrapezoidBounds::operator==(const Acts::SurfaceBounds& sbo) const
{
    // check the type first not to compare apples with oranges
  const Acts::TrapezoidBounds* trabo = dynamic_cast<const Acts::TrapezoidBounds*>(&sbo);
  if (!trabo) return false;
  return (m_boundValues == trabo->m_boundValues);
}

// checking if inside bounds
bool Acts::TrapezoidBounds::inside(const Acts::Vector2D& locpo, double tol1, double tol2) const
{
    if (m_alpha==0.) return insideFull(locpo, tol1, tol2);
    return (insideFull(locpo, tol1, tol2) && !insideExclude(locpo,tol1,tol2));
}

// checking if inside bounds (Full symmetrical Trapezoid)
bool Acts::TrapezoidBounds::insideFull(const Acts::Vector2D& locpo,
                                      double tol1,
                                      double tol2) const
{
  // the cases:
  // the cases:
  double fabsX = fabs(locpo[Acts::eLOC_X]);
  double fabsY = fabs(locpo[Acts::eLOC_Y]);
  // (1) a fast FALSE
  if (fabsY > ( m_boundValues.at(TrapezoidBounds::bv_halfY) + tol2)) return false;
  // (2) a fast FALSE
  if (fabsX > (m_boundValues.at(TrapezoidBounds::bv_maxHalfX) + tol1)) return false;
  // (3) a fast TRUE
  if (fabsX < (m_boundValues.at(TrapezoidBounds::bv_minHalfX) - tol1)) return true;
  // (4) particular case - a rectangle
  if ( m_boundValues.at(TrapezoidBounds::bv_maxHalfX)==m_boundValues.at(TrapezoidBounds::bv_minHalfX) ) return true;
  // (5) /** use basic calculation of a straight line */
  double k = 2.0 * m_boundValues.at(TrapezoidBounds::bv_halfY)/(m_boundValues.at(TrapezoidBounds::bv_maxHalfX) - m_boundValues.at(TrapezoidBounds::bv_minHalfX)) * ( (locpo[Acts::eLOC_X] > 0.) ? 1.0 : -1.0 );
  double d = - fabs(k)*0.5*(m_boundValues.at(TrapezoidBounds::bv_maxHalfX) + m_boundValues.at(TrapezoidBounds::bv_minHalfX));
  return (isAbove(locpo, tol1, tol2, k, d));
}

// checking if local point is inside the exclude area
bool Acts::TrapezoidBounds::insideExclude(const Acts::Vector2D& locpo,
                                         double tol1,
                                         double tol2) const
{

    //line a
    bool alphaBiggerBeta(m_alpha > m_beta);
    double ka   =  alphaBiggerBeta ? tan(M_PI - m_alpha) : tan(m_alpha);
    double kb   =  alphaBiggerBeta ? tan(M_PI - m_beta)  : tan(m_beta);
    double sign =  alphaBiggerBeta ? -1. : 1.;
    double da   =  -m_boundValues.at(TrapezoidBounds::bv_halfY) + sign*ka*m_boundValues.at(TrapezoidBounds::bv_minHalfX);
    double db   =  -m_boundValues.at(TrapezoidBounds::bv_halfY) + sign*kb*m_boundValues.at(TrapezoidBounds::bv_minHalfX);

    return (isAbove(locpo,tol1, tol2,ka,da) && isAbove(locpo,tol1, tol2,kb,db));
}

// checking if local point lies above a line
bool Acts::TrapezoidBounds::isAbove(const Acts::Vector2D& locpo,
                                   double tol1,
                                   double tol2,
                                   double k,
                                   double d) const
{
    // the most tolerant approach for tol1 and tol2
    double sign = k > 0. ?  -1. : + 1.;
    return ( locpo[Acts::eLOC_Y] + tol2 > (k * ( locpo[Acts::eLOC_X] + sign*tol1)+ d) );
}


double Acts::TrapezoidBounds::minDistance(const Acts::Vector2D& pos ) const
{
  const int Np = 4;

  double xl = -m_boundValues.at(TrapezoidBounds::bv_maxHalfX);
  double xr =  m_boundValues.at(TrapezoidBounds::bv_maxHalfX);
  if      (m_alpha !=0.) {
    xl = -m_boundValues.at(TrapezoidBounds::bv_minHalfX)-2.*tan(m_alpha)*m_boundValues.at(TrapezoidBounds::bv_halfY);
  }
  else if (m_beta !=0. ) {
    xr =  m_boundValues.at(TrapezoidBounds::bv_minHalfX)+2.*tan(m_beta )*m_boundValues.at(TrapezoidBounds::bv_halfY);
  }
  double X [4] = { -m_boundValues.at(TrapezoidBounds::bv_minHalfX),
                    xl,
                    xr,
                    m_boundValues.at(TrapezoidBounds::bv_minHalfX)};
  double Y [4] = {  -m_boundValues.at(TrapezoidBounds::bv_halfY),
                    m_boundValues.at(TrapezoidBounds::bv_halfY),
                    m_boundValues.at(TrapezoidBounds::bv_halfY),
                    -m_boundValues.at(TrapezoidBounds::bv_halfY)};

  double dm = 1.e+20;
  double Ao =     0.;
  bool   in =   true;

  for (int i=0; i!=Np; ++i) {

    int j = i+1; if(j==Np) j=0;

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
      if (S <= a ) {
          double d = (A*A)/a;
          if (d<dm) dm=d;
      }
    }
    if (i && in && Ao*A < 0.) in = false;
    Ao = A;
  }
  if (in) return -sqrt(dm);
  else return sqrt(dm);
}

// ostream operator overload
std::ostream& Acts::TrapezoidBounds::dump( std::ostream& sl ) const
{
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(7);
    sl << "Acts::TrapezoidBounds:  (minHlenghtX, maxHlengthX, hlengthY) = " << "("
        << m_boundValues.at(TrapezoidBounds::bv_minHalfX) << ", "
        << m_boundValues.at(TrapezoidBounds::bv_maxHalfX) << ", "
        << m_boundValues.at(TrapezoidBounds::bv_halfY) << ")";
    sl << std::setprecision(-1);
    return sl;
}

