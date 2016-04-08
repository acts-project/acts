///////////////////////////////////////////////////////////////////
// VolumeBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////

#include "Volumes/VolumeBounds.h"

/**Overload of << operator for both, MsgStream and std::ostream for debug output*/ 
MsgStream& Acts::operator << ( MsgStream& sl, const Acts::VolumeBounds& vb)
{ return vb.dump(sl); }

std::ostream& Acts::operator << ( std::ostream& sl, const Acts::VolumeBounds& vb)
{ return vb.dump(sl); }    
