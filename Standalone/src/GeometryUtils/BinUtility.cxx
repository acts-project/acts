
///////////////////////////////////////////////////////////////////
// BinUtility.cxx, ACTS project
///////////////////////////////////////////////////////////////////

//Trk
#include "GeometryUtils/BinUtility.h"
// Gaudi
#include "GaudiKernel/MsgStream.h"
// STD/STL
#include <iostream>

/**Overload of << operator for both, MsgStream and std::ostream for debug output*/ 
MsgStream& Acts::operator << ( MsgStream& sl, const BinUtility& bgen)
{ return bgen.dump(sl); } 

std::ostream& Acts::operator << ( std::ostream& sl, const BinUtility& bgen)
{ return bgen.dump(sl); } 

