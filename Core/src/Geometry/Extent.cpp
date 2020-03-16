// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Extent.hpp"

#include <ostream>

std::ostream& Acts::Extent::toStream(std::ostream& sl) const {
  sl << "Extent in space : " << std::endl;
  for (size_t ib = 0; ib < static_cast<size_t>(binValues); ++ib) {
    sl << "  - value :" << std::setw(10) << binningValueNames[ib]
       << " | range = [" << ranges[ib].first << ", " << ranges[ib].second << "]"
       << std::endl;
  }
  return sl;
}

// Overload of << operator for std::ostream for debug output
std::ostream& Acts::operator<<(std::ostream& sl, const Extent& ext) {
  return ext.toStream(sl);
}