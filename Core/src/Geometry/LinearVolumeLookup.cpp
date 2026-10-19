// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/LinearVolumeLookup.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"

#include <algorithm>

#include "VolumeLookupImpl.hpp"

namespace Acts {
struct LinearVolumeLookup::Impl {
  Config config;
};

LinearVolumeLookup::LinearVolumeLookup() : m_impl(std::make_unique<Impl>()) {}
LinearVolumeLookup::~LinearVolumeLookup() = default;

LinearVolumeLookup::LinearVolumeLookup(const TrackingVolume& world,
                                       const Config& config,
                                       const Logger& logger)
    : m_impl(std::make_unique<Impl>(Impl{config})) {
  const auto warningThreshold = config.warningThreshold;
  if (warningThreshold == 0 || !logger.doPrint(Logging::WARNING)) {
    return;
  }
  std::size_t largest = 0;
  // apply() cannot prune subtrees: explicit child lookups end this scan's
  // scope, and their descendants must not contribute to its warning.
  auto inspect = [&](auto&& self, const TrackingVolume& volume) -> void {
    largest = std::max(
        {largest, volume.volumes().size(), volume.denseVolumes().size()});
    if (const auto confined = volume.confinedVolumes()) {
      for (const auto& child : confined->arrayObjects()) {
        if (child && child->volumeLookup() == nullptr) {
          self(self, *child);
        }
      }
    }
    for (const auto& child : volume.volumes()) {
      if (child.volumeLookup() == nullptr) {
        self(self, child);
      }
    }
  };
  inspect(inspect, world);
  if (largest > warningThreshold) {
    ACTS_WARNING("Linear volume lookup scans up to "
                 << largest << " siblings (warning threshold "
                 << warningThreshold
                 << "). Consider RzKdTree or RzGrid via "
                    "setVolumeLookupFactory on a blueprint node.");
  }
}

const LinearVolumeLookup::Config& LinearVolumeLookup::config() const {
  return m_impl->config;
}

const TrackingVolume* LinearVolumeLookup::lowestTrackingVolume(
    const TrackingVolume& root, const GeometryContext& gctx,
    const Vector3& position, double tolerance) const {
  return detail::linearVolumeScan(root, gctx, position, tolerance);
}
const TrackingVolume* detail::linearVolumeScan(const TrackingVolume& root,
                                               const GeometryContext& gctx,
                                               const Vector3& position,
                                               double tol) {
  if (!root.inside(gctx, position, tol)) {
    return nullptr;
  }

  // confined static volumes - highest hierarchy
  if (const auto confined = root.confinedVolumes()) {
    const TrackingVolume* volume = confined->object(position).get();
    if (volume != nullptr) {
      return volume->lowestTrackingVolume(gctx, position, tol);
    }
  }

  // search for dense volumes
  for (const auto& denseVolume : root.denseVolumes()) {
    if (denseVolume->inside(gctx, position, tol)) {
      return denseVolume.get();
    }
  }

  for (const auto& volume : root.volumes()) {
    if (volume.inside(gctx, position, tol)) {
      return volume.lowestTrackingVolume(gctx, position, tol);
    }
  }

  // there is no lower sub structure
  return &root;
}
}  // namespace Acts
