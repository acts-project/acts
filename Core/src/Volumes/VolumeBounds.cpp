///////////////////////////////////////////////////////////////////
// VolumeBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Volumes/VolumeBounds.hpp"

/**Overload of << operator for std::ostream for debug output*/
std::ostream& Acts::operator << ( std::ostream& sl, const Acts::VolumeBounds& vb)
{ return vb.dump(sl); }
