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
               << " portals=" << m_cfg.portals << " sensitives="
               << m_cfg.sensitives << " passives=" << m_cfg.passives);
}

TryAllNavigationPolicy::TryAllNavigationPolicy(const GeometryContext& gctx,
                                               const TrackingVolume& volume,
                                               const Logger& logger)
    : TryAllNavigationPolicy(gctx, volume, logger, {}) {}

void TryAllNavigationPolicy::initializeCandidates(
    [[maybe_unused]] const GeometryContext& gctx,
    const NavigationArguments& args, NavigationPolicyState& /*state*/,
    AppendOnlyNavigationStream& stream, const Logger& logger) const {
  ACTS_VERBOSE("TryAllNavigationPolicy initializing candidates for volume "
               << m_volume->volumeName());
  ACTS_VERBOSE("~> Config: portals=" << m_cfg.portals
                                     << " sensitives=" << m_cfg.sensitives
                                     << " passives=" << m_cfg.passives);
  assert(m_volume != nullptr);

  std::size_t numCandidates = 0;

  if (m_cfg.portals) {
    for (const auto& portal : m_volume->portals()) {
      stream.addPortalCandidate(portal);
      numCandidates++;
    }
  }

  if (!(m_cfg.sensitives || m_cfg.passives)) {
    return;
  }

  for (const auto& surface : m_volume->surfaces()) {
    bool isSensitive = surface.isSensitive();
    if ((m_cfg.passives && !isSensitive) || (m_cfg.sensitives && isSensitive)) {
      stream.addSurfaceCandidate(surface, args.tolerance);
      numCandidates++;
    }
  }

  ACTS_VERBOSE("TryAllNavigationPolicy added " << numCandidates
                                               << " candidates to the stream");
}

void TryAllNavigationPolicy::connect(NavigationDelegate& delegate) const {
  connectDefault<TryAllNavigationPolicy>(delegate);
}

}  // namespace Acts
