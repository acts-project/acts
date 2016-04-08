///////////////////////////////////////////////////////////////////
// GeometryIDCalculator.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "GeometryUtils/GeometryID.h"

bool Acts::operator< ( const Acts::GeometryID& one, const Acts::GeometryID& two )
{ return (one.value() < two.value());  }

bool Acts::operator<=( const Acts::GeometryID& one, const Acts::GeometryID& two )
{ return (one.value() <= two.value()); }

bool Acts::operator> ( const Acts::GeometryID& one, const Acts::GeometryID& two )
{ return (one.value() > two.value());  } 

bool Acts::operator>=( const Acts::GeometryID& one, const Acts::GeometryID& two )
{ return (one.value() >= two.value()); }

MsgStream& Acts::operator << ( MsgStream& sl, const Acts::GeometryID& tid)
{
    sl << tid.value();
    return sl;
}

std::ostream& Acts::operator << ( std::ostream& sl, const Acts::GeometryID& tid)
{
    sl << tid.value();
    return sl;
}
