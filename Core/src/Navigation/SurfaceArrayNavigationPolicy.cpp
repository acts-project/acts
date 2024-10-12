// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"

namespace Acts {

SurfaceArrayNavigationPolicy::SurfaceArrayNavigationPolicy(
    const GeometryContext& gctx, const TrackingVolume& volume,
    const Logger& logger, Config config) {
  ACTS_VERBOSE("Constructing SurfaceArrayNavigationPolicy for volume "
               << volume.volumeName());
  ACTS_VERBOSE("Layer type is " << config.layerType);
  // collect sensitive surfaces from the volume
  // std::vector<std::shared_ptr<const Surface>>;
}

void SurfaceArrayNavigationPolicy::updateState(
    const NavigationArguments& args) const {}

void SurfaceArrayNavigationPolicy::connect(NavigationDelegate& delegate) const {
  connectDefault<SurfaceArrayNavigationPolicy>(delegate);
}

}  // namespace Acts
