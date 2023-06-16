// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Direction.hpp"

#include <cstdlib>

std::string Acts::Direction::toString() const {
  switch (m_value) {
    case Value::Positive:
      return "forward";
    case Value::Negative:
      return "backward";
    default:
      assert(false && "Invalid direction (shouldn't happen)");
      std::abort();
  }
}

std::ostream& Acts::operator<<(std::ostream& os, Direction dir) {
  os << dir.toString();
  return os;
}
