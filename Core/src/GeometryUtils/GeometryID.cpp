///////////////////////////////////////////////////////////////////
// GeometryIDCalculator.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Utilities/GeometryID.hpp"

bool Acts::operator< ( const Acts::GeometryID& one, const Acts::GeometryID& two )
{ return (one.value() < two.value());  }

bool Acts::operator<=( const Acts::GeometryID& one, const Acts::GeometryID& two )
{ return (one.value() <= two.value()); }

bool Acts::operator> ( const Acts::GeometryID& one, const Acts::GeometryID& two )
{ return (one.value() > two.value());  }

bool Acts::operator>=( const Acts::GeometryID& one, const Acts::GeometryID& two )
{ return (one.value() >= two.value()); }

std::ostream& Acts::operator << ( std::ostream& sl, const Acts::GeometryID& tid)
{
    sl << tid.value();
    return sl;
}
