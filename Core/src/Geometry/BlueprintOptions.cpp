// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/BlueprintOptions.hpp"

#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/CylinderNavigationPolicy.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"

#include <memory>

namespace Acts::Experimental {

void BlueprintOptions::validate() const {
  if (!defaultNavigationPolicyFactory) {
    throw std::invalid_argument("Navigation policy factory is nullptr");
  }
}

std::unique_ptr<NavigationPolicyFactory>
BlueprintOptions::makeDefaultNavigationPolicyFactory() {
  return NavigationPolicyFactory{}
      .add([](const GeometryContext& gctx, const TrackingVolume& volume,
              const Logger& logger) -> std::unique_ptr<INavigationPolicy> {
        // The cylinder policy only produces portal candidates, so it can only
        // be the sole policy of a volume that does not confine any surfaces.
        if (volume.surfaces().empty() &&
            CylinderNavigationPolicy::isApplicable(volume, logger)) {
          return std::make_unique<CylinderNavigationPolicy>(gctx, volume,
                                                            logger);
        }
        // Fall back to try-all for every volume the cylinder policy cannot
        // describe on its own.
        return std::make_unique<TryAllNavigationPolicy>(gctx, volume, logger);
      })
      .asUniquePtr();
}

}  // namespace Acts::Experimental
