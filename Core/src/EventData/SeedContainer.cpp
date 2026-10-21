// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SeedContainer.hpp"

#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <limits>

namespace Acts {

static_assert(std::ranges::random_access_range<SeedContainer>);

SeedContainer::SeedContainer() noexcept = default;

SeedContainer::SeedContainer(const SeedContainer &other) noexcept = default;

SeedContainer::SeedContainer(SeedContainer &&other) noexcept = default;

SeedContainer &SeedContainer::operator=(const SeedContainer &other) noexcept =
    default;

SeedContainer &SeedContainer::operator=(SeedContainer &&other) noexcept =
    default;

void SeedContainer::reserve(Index size, float averageSpacePoints) noexcept {
  m_spacePointOffsets.reserve(size);
  m_spacePointCounts.reserve(size);
  m_qualities.reserve(size);
  m_vertexZs.reserve(size);
  m_spacePoints.reserve(static_cast<std::size_t>(size * averageSpacePoints));
}

void SeedContainer::clear() noexcept {
  m_size = 0;

  m_spacePointOffsets.clear();
  m_spacePointCounts.clear();
  m_qualities.clear();
  m_vertexZs.clear();
  m_spacePoints.clear();
}

void SeedContainer::assignSpacePointContainer(
    SpacePointContainer &&spacePointContainer) noexcept {
  auto movedContainer =
      std::make_shared<SpacePointContainer>(std::move(spacePointContainer));

  m_sharedConstSpacePointContainer = movedContainer;
  m_constSpacePointContainer = movedContainer.get();
  m_mutableSpacePointContainer = movedContainer.get();
}

void SeedContainer::assignSpacePointContainer(
    SpacePointContainer &spacePointContainer) noexcept {
  m_sharedConstSpacePointContainer = nullptr;
  m_mutableSpacePointContainer = &spacePointContainer;
  m_constSpacePointContainer = &spacePointContainer;
}

void SeedContainer::assignSpacePointContainer(
    const SpacePointContainer &spacePointContainer) noexcept {
  m_sharedConstSpacePointContainer = nullptr;
  m_constSpacePointContainer = &spacePointContainer;
  m_mutableSpacePointContainer = nullptr;
}

void SeedContainer::assignSpacePointContainer(
    const std::shared_ptr<SpacePointContainer> &spacePointContainer) noexcept {
  m_sharedConstSpacePointContainer = spacePointContainer;
  m_constSpacePointContainer = spacePointContainer.get();
  m_mutableSpacePointContainer = spacePointContainer.get();
}

void SeedContainer::assignSpacePointContainer(
    const std::shared_ptr<const SpacePointContainer>
        &spacePointContainer) noexcept {
  m_sharedConstSpacePointContainer = spacePointContainer;
  m_constSpacePointContainer = spacePointContainer.get();
  m_mutableSpacePointContainer = nullptr;
}

bool SeedContainer::hasSpacePointContainer() const noexcept {
  return m_constSpacePointContainer != nullptr;
}

bool SeedContainer::hasMutableSpacePointContainer() const noexcept {
  return m_mutableSpacePointContainer != nullptr;
}

SpacePointContainer &SeedContainer::mutableSpacePointContainer() {
  if (!hasMutableSpacePointContainer()) {
    throw std::logic_error(
        "No mutable SpacePointContainer assigned to SeedContainer");
  }
  return *m_mutableSpacePointContainer;
}

const SpacePointContainer &SeedContainer::spacePointContainer() const {
  if (!hasSpacePointContainer()) {
    throw std::logic_error("No SpacePointContainer assigned to SeedContainer");
  }
  return *m_constSpacePointContainer;
}

MutableSeedProxy SeedContainer::createSeed() noexcept {
  ++m_size;

  m_spacePointOffsets.push_back(
      static_cast<std::uint32_t>(m_spacePoints.size()));
  m_spacePointCounts.push_back(static_cast<std::uint8_t>(0));
  m_qualities.push_back(-std::numeric_limits<float>::infinity());
  m_vertexZs.push_back(0.f);

  return MutableProxy(*this, size() - 1);
}

void SeedContainer::copyFrom(Index index, const SeedContainer &otherContainer,
                             Index otherIndex, SeedColumns columnsToCopy) {
  if (index >= size()) {
    throw std::out_of_range("Index out of range in SeedContainer");
  }
  if (otherIndex >= otherContainer.size()) {
    throw std::out_of_range("Other index out of range in SeedContainer");
  }

  if (ACTS_CHECK_BIT(columnsToCopy, SeedColumns::SpacePointIndices)) {
    at(index).assignSpacePointIndices(
        otherContainer.at(otherIndex).spacePointIndices());
  }
  if (ACTS_CHECK_BIT(columnsToCopy, SeedColumns::Quality)) {
    at(index).quality() = otherContainer.at(otherIndex).quality();
  }
  if (ACTS_CHECK_BIT(columnsToCopy, SeedColumns::VertexZ)) {
    at(index).vertexZ() = otherContainer.at(otherIndex).vertexZ();
  }
}

}  // namespace Acts
