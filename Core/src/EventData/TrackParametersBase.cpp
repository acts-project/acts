// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
