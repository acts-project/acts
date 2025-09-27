// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SeedContainer2.hpp"

#include "Acts/EventData/Types.hpp"

namespace Acts {

SeedContainer2::SeedContainer2() noexcept = default;

SeedContainer2::SeedContainer2(const SeedContainer2 &other) noexcept
    : m_size(other.m_size) {
  knownColumns() = other.knownColumns();
}

SeedContainer2::SeedContainer2(SeedContainer2 &&other) noexcept
    : m_size(other.m_size) {
  knownColumns() = std::move(other).knownColumns();

  other.m_size = 0;
}

SeedContainer2 &SeedContainer2::operator=(
    const SeedContainer2 &other) noexcept {
  if (this == &other) {
    return *this;
  }

  m_size = other.m_size;
  knownColumns() = other.knownColumns();

  return *this;
}

SeedContainer2 &SeedContainer2::operator=(SeedContainer2 &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  m_size = other.m_size;
  knownColumns() = std::move(other).knownColumns();

  other.m_size = 0;

  return *this;
}

void SeedContainer2::reserve(std::size_t size,
                             float averageSpacePoints) noexcept {
  m_spacePointOffsets.reserve(size);
  m_spacePointCounts.reserve(size);
  m_qualities.reserve(size);
  m_vertexZs.reserve(size);
  m_spacePoints.reserve(static_cast<std::size_t>(size * averageSpacePoints));
}

void SeedContainer2::clear() noexcept {
  m_size = 0;

  m_spacePointOffsets.clear();
  m_spacePointCounts.clear();
  m_qualities.clear();
  m_vertexZs.clear();
  m_spacePoints.clear();
}

void SeedContainer2::assignSpacePointIndices(
    Index index, std::span<const SpacePointIndex2> spacePointIndices) {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SpacePointContainer2");
  }
  if (m_spacePointCounts[index] != 0) {
    throw std::logic_error("Space points already assigned to the seed");
  }

  m_spacePointOffsets[index] =
      static_cast<SpacePointIndex2>(m_spacePoints.size());
  m_spacePointCounts[index] =
      static_cast<std::uint8_t>(spacePointIndices.size());
  m_spacePoints.insert(m_spacePoints.end(), spacePointIndices.begin(),
                       spacePointIndices.end());
}

}  // namespace Acts
