///////////////////////////////////////////////////////////////////
// SurfaceBounds.cxx, ACTS project
///////////////////////////////////////////////////////////////////


//Trk
#include "Surfaces/SurfaceBounds.h"

/**Overload of << operator for both, MsgStream and std::ostream for debug output*/ 
MsgStream& Acts::operator << ( MsgStream& sl, const SurfaceBounds& sb)
{ return sb.dump(sl); }

std::ostream& Acts::operator << ( std::ostream& sl, const SurfaceBounds& sb)
{ return sb.dump(sl); }    
