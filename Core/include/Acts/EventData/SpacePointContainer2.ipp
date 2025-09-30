// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"

#include "Acts/EventData/SpacePointProxy2.hpp"

namespace Acts {

inline MutableSpacePointProxy2 SpacePointContainer2::at(Index index) {
  if (index >= size()) {
    throw std::out_of_range(
        "Index out of range in SpacePointContainer2: " + std::to_string(index) +
        " >= " + std::to_string(size()));
  }
  return MutableProxy(*this, index);
}

inline ConstSpacePointProxy2 SpacePointContainer2::at(Index index) const {
  if (index >= size()) {
    throw std::out_of_range(
        "Index out of range in SpacePointContainer2: " + std::to_string(index) +
        " >= " + std::to_string(size()));
  }
  return ConstProxy(*this, index);
}

inline MutableSpacePointProxy2 SpacePointContainer2::operator[](
    Index index) noexcept {
  return MutableProxy(*this, index);
}

inline ConstSpacePointProxy2 SpacePointContainer2::operator[](
    Index index) const noexcept {
  return ConstProxy(*this, index);
}

}  // namespace Acts
