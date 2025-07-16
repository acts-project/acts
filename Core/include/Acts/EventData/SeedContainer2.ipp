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

namespace Acts::Experimental {

inline MutableSeedProxy2 SeedContainer2::createSeed(
    std::span<const SpacePointIndex2> spacePoints) noexcept {
  m_sizes.push_back(static_cast<std::uint8_t>(spacePoints.size()));
  m_spacePointOffsets.push_back(m_spacePoints.size());
  m_qualities.push_back(-std::numeric_limits<float>::infinity());
  m_vertexZs.push_back(0.f);
  m_spacePoints.insert(m_spacePoints.end(), spacePoints.begin(),
                       spacePoints.end());
  return MutableProxyType(*this, size() - 1);
}

inline MutableSeedProxy2 SeedContainer2::at(IndexType index) {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SeedContainer2");
  }
  return MutableProxyType(*this, index);
}

inline ConstSeedProxy2 SeedContainer2::at(IndexType index) const {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SeedContainer2");
  }
  return ConstProxyType(*this, index);
}

}  // namespace Acts::Experimental
