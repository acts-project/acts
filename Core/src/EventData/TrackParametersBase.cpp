// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// STL include(s)
#include <iomanip>

// Acts include(s)
#include "Acts/EventData/TrackParametersBase.hpp"

namespace Acts {
std::ostream&
TrackParametersBase::print(std::ostream& sl) const
{
  // set stream output format
  auto old_precision = sl.precision(7);
  auto old_flags     = sl.setf(std::ios::fixed);

  sl << " * TrackParameters:" << std::endl;
  sl << parameters() << std::endl;
  sl << " * charge: " << charge() << std::endl;
  if (covariance() != nullptr) {
    sl << " * covariance matrix = " << *covariance() << std::endl;
  } else {
    sl << " * covariance matrix = " << covariance() << std::endl;
  }
  sl << " * corresponding global parameters:" << std::endl;
  sl << " *    position  (x,  y,  z ) = (" << position().x() << ", "
     << position().y() << ", " << position().z() << ")" << std::endl;
  sl << " *    momentum  (px, py, pz) = (" << momentum().x() << ", "
     << momentum().y() << ", " << momentum().z() << ")" << std::endl;

  // reset stream format
  sl.precision(old_precision);
  sl.setf(old_flags);

  return sl;
}
}  // namespace Acts
