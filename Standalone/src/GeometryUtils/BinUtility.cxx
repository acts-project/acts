
///////////////////////////////////////////////////////////////////
// BinUtility.cxx, ACTS project
///////////////////////////////////////////////////////////////////

//Trk
#include "ACTS/Utilities/BinUtility.h"

// STD/STL
#include <iostream>

/**Overload of << operator for std::ostream for debug output*/
std::ostream& Acts::operator << ( std::ostream& sl, const BinUtility& bgen)
{ return bgen.dump(sl); }

