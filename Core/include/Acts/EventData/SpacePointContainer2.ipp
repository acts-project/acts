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

namespace Acts::Experimental {

inline MutableSpacePointProxy2
SpacePointContainer2::createSpacePoint() noexcept {
  for (const auto &[name, column] : m_namedColumns) {
    column.first->emplace_back();
  }

  return MutableProxyType(*this, size() - 1);
}

inline MutableSpacePointProxy2 SpacePointContainer2::at(IndexType index) {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SpacePointContainer2");
  }
  return MutableProxyType(*this, index);
}

inline ConstSpacePointProxy2 SpacePointContainer2::at(IndexType index) const {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SpacePointContainer2");
  }
  return ConstProxyType(*this, index);
}

inline MutableSpacePointProxy2 SpacePointContainer2::operator[](
    IndexType index) noexcept {
  return MutableProxyType(*this, index);
}

inline ConstSpacePointProxy2 SpacePointContainer2::operator[](
    IndexType index) const noexcept {
  return ConstProxyType(*this, index);
}

}  // namespace Acts::Experimental
