// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/TryAllNavigationPolicies.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/NavigationStream.hpp"

namespace Acts {

TryAllPortalNavigationPolicy::TryAllPortalNavigationPolicy(
    const GeometryContext& /*gctx*/, const TrackingVolume& volume,
    const Logger& /*logger*/)
    : m_volume(&volume) {}

void TryAllPortalNavigationPolicy::updateState(
    const NavigationArguments& args) const {
  assert(m_volume != nullptr);

  for (const auto& portal : m_volume->portals()) {
    args.main.addPortalCandidate(portal);
  };
}

void TryAllPortalNavigationPolicy::connect(NavigationDelegate& delegate) const {
  connectDefault<TryAllPortalNavigationPolicy>(delegate);
}

TryAllSurfaceNavigationPolicy::TryAllSurfaceNavigationPolicy(
    const GeometryContext& /*gctx*/, const TrackingVolume& volume,
    const Logger& /*logger*/)
    : m_volume(&volume) {}

void TryAllSurfaceNavigationPolicy::updateState(
    const NavigationArguments& args) const {
  assert(m_volume != nullptr);

  for (const auto& surface : m_volume->surfaces()) {
    args.main.addSurfaceCandidate(surface, args.tolerance);
  };
}

void TryAllSurfaceNavigationPolicy::connect(
    NavigationDelegate& delegate) const {
  connectDefault<TryAllSurfaceNavigationPolicy>(delegate);
}

}  // namespace Acts
