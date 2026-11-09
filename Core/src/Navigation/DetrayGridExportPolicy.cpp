// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/DetrayGridExportPolicy.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace Acts {

namespace {

// Composite policies only visit their direct children, so recurse explicitly
void collectSurfaceArrayPolicies(
    const INavigationPolicy& policy,
    std::vector<const SurfaceArrayNavigationPolicy*>& found) {
  policy.visit([&](const INavigationPolicy& visited) {
    if (&visited != &policy) {
      collectSurfaceArrayPolicies(visited, found);
      return;
    }
    if (const auto* surfaceArray =
            dynamic_cast<const SurfaceArrayNavigationPolicy*>(&visited);
        surfaceArray != nullptr) {
      found.push_back(surfaceArray);
    }
  });
}

}  // namespace

DetrayGridExportPolicy::DetrayGridExportPolicy(
    const GeometryContext& /*gctx*/, const TrackingVolume& volume,
    const Logger& logger, std::span<const INavigationPolicy* const> siblings,
    const Config& config)
    : m_cfg{config} {
  std::vector<const SurfaceArrayNavigationPolicy*> found;
  for (const INavigationPolicy* sibling : siblings) {
    if (sibling != nullptr) {
      collectSurfaceArrayPolicies(*sibling, found);
    }
  }

  if (found.size() != 1) {
    throw std::invalid_argument(
        "DetrayGridExportPolicy in volume '" + volume.volumeName() +
        "' needs exactly one SurfaceArrayNavigationPolicy added before it to "
        "the same navigation policy factory, found " +
        std::to_string(found.size()));
  }
  m_surfaceArrayPolicy = found.front();

  ACTS_VERBOSE("DetrayGridExportPolicy created for volume "
               << volume.volumeName()
               << " with vertexZHalfRange=" << m_cfg.vertexZHalfRange
               << " minPt=" << m_cfg.minPt << " bField=" << m_cfg.bField
               << " fixedWindow=" << m_cfg.fixedWindow);
}

void DetrayGridExportPolicy::initializeCandidates(
    const GeometryContext& /*gctx*/, const NavigationArguments& /*args*/,
    NavigationPolicyState& /*state*/, AppendOnlyNavigationStream& /*stream*/,
    const Logger& /*logger*/) const {
  // Export information only: the sibling provides the candidates
}

void DetrayGridExportPolicy::connect(NavigationDelegate& delegate) const {
  delegate.connect<&INavigationPolicy::noopInitializeCandidates>();
}

const DetrayGridExportPolicy::Config& DetrayGridExportPolicy::config() const {
  return m_cfg;
}

const SurfaceArrayNavigationPolicy& DetrayGridExportPolicy::surfaceArrayPolicy()
    const {
  return *m_surfaceArrayPolicy;
}

}  // namespace Acts
