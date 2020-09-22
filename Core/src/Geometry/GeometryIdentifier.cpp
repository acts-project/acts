// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <iomanip>
#include <ostream>

std::ostream& Acts::operator<<(std::ostream& os, Acts::GeometryIdentifier id) {
  os << "[ " << std::setw(3) << id.volume();
  os << " | " << std::setw(3) << id.boundary();
  os << " | " << std::setw(3) << id.layer();
  os << " | " << std::setw(3) << id.approach();
  os << " | " << std::setw(4) << id.sensitive() << " ]";
  return os;
}
