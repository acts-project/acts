///////////////////////////////////////////////////////////////////
// TriangleBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "Surfaces/TriangleBounds.h"
// STD/STL
#include <iostream>
#include <iomanip>

// default constructor
Acts::TriangleBounds::TriangleBounds() :
    Acts::PlanarBounds(),
    m_boundValues(TriangleBounds::bv_length,0.)
{}

// rectangle constructor - float constructor
Acts::TriangleBounds::TriangleBounds(std::vector<std::pair<float,float> > vertices) :
    Acts::PlanarBounds(),
    m_boundValues(TriangleBounds::bv_length,0.)
{
   std::vector< std::pair<float, float> >::iterator vIt  = vertices.begin();
   std::vector< std::pair<float, float> >::iterator vItE = vertices.end();
   for (size_t ib = 0 ; vIt != vItE; ++vIt, ++ib){
       m_boundValues[2*ib]    = (*vIt).first;
       m_boundValues[2*ib+1]  = (*vIt).second;
       if (ib==2) break;
   }
}

// rectangle constructor - double constructor
Acts::TriangleBounds::TriangleBounds(std::vector<std::pair<double,double> > vertices) :
    Acts::PlanarBounds(),
    m_boundValues(TriangleBounds::bv_length,0.)
{
    std::vector< std::pair<double, double> >::iterator vIt  = vertices.begin();
    std::vector< std::pair<double, double> >::iterator vItE = vertices.end();
    for (size_t ib = 0 ; vIt != vItE; ++vIt, ++ib){
        m_boundValues[2*ib]    = (*vIt).first;
        m_boundValues[2*ib+1]  = (*vIt).second;
        if (ib==2) break;
    }
}

// constructor from three points
Acts::TriangleBounds::TriangleBounds( const Acts::Vector2D& p1, const Acts::Vector2D& p2, const Acts::Vector2D& p3) :
    Acts::PlanarBounds(),
    m_boundValues(TriangleBounds::bv_length,0.)
{
    m_boundValues[TriangleBounds::bv_x1] = p1.x();
    m_boundValues[TriangleBounds::bv_y1] = p1.y();
    m_boundValues[TriangleBounds::bv_x2] = p2.x();
    m_boundValues[TriangleBounds::bv_y2] = p2.y();
    m_boundValues[TriangleBounds::bv_x3] = p3.x();
    m_boundValues[TriangleBounds::bv_y3] = p3.y();
}

// copy constructor
Acts::TriangleBounds::TriangleBounds(const TriangleBounds& tribo) :
  Acts::PlanarBounds(),
  m_boundValues(tribo.m_boundValues)
{}

// destructor
Acts::TriangleBounds::~TriangleBounds()
{}

Acts::TriangleBounds& Acts::TriangleBounds::operator=(const TriangleBounds& tribo)
{
  if (this!=&tribo)
      m_boundValues = tribo.m_boundValues;
  return *this;
}

bool Acts::TriangleBounds::operator==(const Acts::SurfaceBounds& sbo) const
{
  // check the type first not to compare apples with oranges
  const Acts::TriangleBounds* tribo = dynamic_cast<const Acts::TriangleBounds*>(&sbo);
  if (!tribo) return false;
  return (m_boundValues == tribo->m_boundValues);
}

double Acts::TriangleBounds::minDistance(const Acts::Vector2D& pos ) const
{
  const int Np = 3;

  double X [3] = { m_boundValues[TriangleBounds::bv_x1] ,
                   m_boundValues[TriangleBounds::bv_x2] ,
                   m_boundValues[TriangleBounds::bv_x3] };
  double Y [3] = { m_boundValues[TriangleBounds::bv_y1] ,
                   m_boundValues[TriangleBounds::bv_y2] ,
                   m_boundValues[TriangleBounds::bv_y3] };

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
std::ostream& Acts::TriangleBounds::dump( std::ostream& sl ) const
{
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(7);
    sl << "Acts::TriangleBounds:  generating vertices (X, Y)";
    sl <<  "(" << m_boundValues[TriangleBounds::bv_x1] << " , " << m_boundValues[TriangleBounds::bv_y1] << ") " << '\n';
    sl <<  "(" << m_boundValues[TriangleBounds::bv_x2] << " , " << m_boundValues[TriangleBounds::bv_y2] << ") " << '\n';
    sl <<  "(" << m_boundValues[TriangleBounds::bv_x3] << " , " << m_boundValues[TriangleBounds::bv_y3] << ") ";
    sl << std::setprecision(-1);
    return sl;
}

