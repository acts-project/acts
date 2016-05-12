// STL include(s)
#include <iostream>
#include <iomanip>

// ACTS include(s)
#include "ACTS/EventData/TrackParametersBase.hpp"

namespace Acts
{
  std::ostream& TrackParametersBase::dump(std::ostream& sl) const
  {
    sl << std::setiosflags(std::ios::fixed);
    sl << std::setprecision(7);
    sl << " * TrackParameters:" << std::endl;
    sl << parameters() << std::endl;
    sl << " * charge: " << charge() << std::endl;
    sl << " * covariance matrix = " << covariance() << std::endl;
    sl << " * corresponding global parameters:" << std::endl;
    sl << " *    position  (x,  y,  z ) = ("
       << position().x() << ", "
       << position().y() << ", "
       << position().z() << ")" << std::endl;
    sl << " *    momentum  (px, py, pz) = ("
       << momentum().x() << ", "
       << momentum().y() << ", "
       << momentum().z() << ")" << std::endl;
    sl << std::setprecision(-1);

    return sl;
  }

  std::ostream& operator<<(std::ostream& sl,const TrackParametersBase& p)
  {
    return p.dump(sl);
  }
} // end of namespace Acts
