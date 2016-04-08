///////////////////////////////////////////////////////////////////
// VolumeBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "Volumes/VolumeBounds.h"

/**Overload of << operator for std::ostream for debug output*/
std::ostream& Acts::operator << ( std::ostream& sl, const Acts::VolumeBounds& vb)
{ return vb.dump(sl); }
