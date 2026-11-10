// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Detray/DetraySurfaceGridPolicy.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"

namespace ActsPlugins {

DetraySurfaceGridPolicy::DetraySurfaceGridPolicy(
    const Acts::GeometryContext& /*gctx*/,
    const Acts::TrackingVolume& /*volume*/, const Acts::Logger& /*logger*/,
    const DetraySurfaceGridConfig& config)
    : m_cfg(config) {}

DetraySurfaceGridPolicy::DetraySurfaceGridPolicy(
    const Acts::GeometryContext& gctx, const Acts::TrackingVolume& volume,
    const Acts::Logger& logger)
    : DetraySurfaceGridPolicy(gctx, volume, logger, DetraySurfaceGridConfig{}) {
}

void DetraySurfaceGridPolicy::initializeCandidates(
    const Acts::GeometryContext& /*gctx*/,
    const Acts::NavigationArguments& /*args*/,
    Acts::NavigationPolicyState& /*state*/,
    Acts::AppendOnlyNavigationStream& /*stream*/,
    const Acts::Logger& /*logger*/) const {}

void DetraySurfaceGridPolicy::connect(
    Acts::NavigationDelegate& delegate) const {
  connectDefault<DetraySurfaceGridPolicy>(delegate);
}

const DetraySurfaceGridConfig& DetraySurfaceGridPolicy::config() const {
  return m_cfg;
}

}  // namespace ActsPlugins
