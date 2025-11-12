// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/TryAllNavigationPolicy.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/NavigationStream.hpp"

namespace Acts {

TryAllNavigationPolicy::TryAllNavigationPolicy(const GeometryContext& /*gctx*/,
                                               const TrackingVolume& volume,
                                               const Logger& logger,
                                               const Config& config)
    : m_cfg{config}, m_volume(&volume) {
  assert(m_volume != nullptr);
  ACTS_VERBOSE("TryAllNavigationPolicy created for volume "
               << m_volume->volumeName() << " with config: "
               << " portals=" << m_cfg.resolvePortales
               << " sensitives=" << m_cfg.resolveSensitives
               << " passives=" << m_cfg.resolvePassives);
}

TryAllNavigationPolicy::TryAllNavigationPolicy(const GeometryContext& gctx,
                                               const TrackingVolume& volume,
                                               const Logger& logger)
    : TryAllNavigationPolicy(gctx, volume, logger, {}) {}

void TryAllNavigationPolicy::initializeCandidates(
    const NavigationArguments& args, AppendOnlyNavigationStream& stream,
    const Logger& logger) const {
  ACTS_VERBOSE("TryAllNavigationPolicy");
  assert(m_volume != nullptr);

  if (m_cfg.portals) {
    for (const auto& portal : m_volume->portals()) {
      stream.addPortalCandidate(portal);
    }
  }

  if (!(m_cfg.resolveSensitives || m_cfg.resolvePassives ||
        m_cfg.resolveMaterial)) {
    return;
  }

  for (const auto& surface : m_volume->surfaces()) {
    bool isSensitive = surface.associatedDetectorElement() != nullptr;
    bool hasMaterial = surface.surfaceMaterial() != nullptr;
    if ((m_cfg.resolvePassives && !isSensitive) ||
        (m_cfg.resolveSensitives && isSensitive) ||
        (m_cfg.resolveMaterial && hasMaterial)) {
      stream.addSurfaceCandidate(surface, args.tolerance);
    }
  }
}

void TryAllNavigationPolicy::connect(NavigationDelegate& delegate) const {
  connectDefault<TryAllNavigationPolicy>(delegate);
}

}  // namespace Acts
