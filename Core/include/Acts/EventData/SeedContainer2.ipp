// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SeedContainer2.hpp"

#include "Acts/EventData/SeedProxy2.hpp"

#include <limits>

namespace Acts {

inline MutableSeedProxy2 SeedContainer2::createSeed() noexcept {
  ++m_size;

  m_spacePointOffsets.push_back(m_spacePoints.size());
  m_spacePointCounts.push_back(static_cast<std::uint8_t>(0));
  m_qualities.push_back(-std::numeric_limits<float>::infinity());
  m_vertexZs.push_back(0.f);

  return MutableProxy(*this, size() - 1);
}

inline MutableSeedProxy2 SeedContainer2::at(Index index) {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SeedContainer2");
  }
  return MutableProxy(*this, index);
}

inline ConstSeedProxy2 SeedContainer2::at(Index index) const {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SeedContainer2");
  }
  return ConstProxy(*this, index);
}

inline MutableSeedProxy2 SeedContainer2::operator[](Index index) noexcept {
  return MutableProxy(*this, index);
}

inline ConstSeedProxy2 SeedContainer2::operator[](Index index) const noexcept {
  return ConstProxy(*this, index);
}

}  // namespace Acts
