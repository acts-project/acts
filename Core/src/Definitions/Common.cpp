// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Common.hpp"

#include <cassert>
#include <cstdlib>
#include <ostream>

std::ostream& Acts::operator<<(std::ostream& os, MaterialUpdateMode mode) {
  switch (mode) {
    case MaterialUpdateMode::NoUpdate:
      os << "None (0)";
      break;
    case MaterialUpdateMode::PreUpdate:
      os << "PreUpdate (1)";
      break;
    case MaterialUpdateMode::PostUpdate:
      os << "PostUpdate (2)";
      break;
    case MaterialUpdateMode::FullUpdate:
      os << "FullUpdate (3)";
      break;
    default:
      assert(false && "Invalid material update mode (shouldn't happen)");
      std::abort();
  }
  return os;
}
