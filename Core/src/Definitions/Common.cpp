// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Common.hpp"

#include <ostream>

std::string Acts::toString(Direction dir) {
  switch (dir) {
    case Direction::Forward:
      return "forward";
    case Direction::Backward:
      return "backward";
    default:
      assert(false && "Invalid direction (shouldn't happen)");
      std::abort();
  }
}

std::ostream& Acts::operator<<(std::ostream& os, Direction dir) {
  os << toString(dir);
  return os;
}

std::ostream& Acts::operator<<(std::ostream& os,
                               MaterialUpdateStage matUpdate) {
  switch (matUpdate) {
    case MaterialUpdateStage::PreUpdate:
      os << "PreUpdate (-1)";
      break;
    case MaterialUpdateStage::PostUpdate:
      os << "PostUpdate (1)";
      break;
    case MaterialUpdateStage::FullUpdate:
      os << "FullUpdate (0)";
      break;
    default:
      assert(false && "Invalid material update stage (shouldn't happen)");
      std::abort();
  }
  return os;
}
