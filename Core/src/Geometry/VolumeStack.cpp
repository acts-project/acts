// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/VolumeStack.hpp"

namespace Acts {

VolumeStack::VolumeStack(
    std::vector<Volume*>& volumes, AxisDirection direction,
    std::pair<VolumeResizeStrategy, VolumeResizeStrategy> resizeStrategies)
    : Volume(initialVolume(volumes)),
      m_direction(direction),
      m_resizeStrategies{resizeStrategies.first, resizeStrategies.second},
      m_volumes(volumes) {}

Volume& VolumeStack::initialVolume(std::span<Volume*> volumes) {
  if (volumes.empty()) {
    throw std::invalid_argument("VolumeStack requires at least one volume");
  }
  return *volumes.front();
}

bool VolumeStack::isGapVolume(const Volume& volume) const {
  return std::ranges::any_of(
      gaps(), [&](const auto& gap) { return gap.get() == &volume; });
}

std::shared_ptr<Volume> VolumeStack::addGapVolume(
    const Transform3& transform, const std::shared_ptr<VolumeBounds>& bounds) {
  auto gapVolume = std::make_shared<Volume>(transform, bounds);
  m_gaps.push_back(gapVolume);
  return gapVolume;
}

const std::vector<std::shared_ptr<Volume>>& VolumeStack::gaps() const {
  return m_gaps;
}

}  // namespace Acts
