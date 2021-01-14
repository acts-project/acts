// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/PlanarSurfaceDrift.hpp"

ActsFatras::PlanarSurfaceDrift::Segment2D
ActsFatras::PlanarSurfaceDrift::toReadout(const Acts::GeometryContext& gctx,
                                          const Acts::Surface& surface,
                                          double depletion,
                                          const Acts::Vector3& pos,
                                          const Acts::Vector3& dir,
                                          const Acts::Vector3& driftDir) const {
  // Transform the hit & direction into the local surface frame
  const auto& invTransform = surface.transform(gctx).inverse();
  Acts::Vector2 pos2Local = (invTransform * pos).segment<2>(0);
  Acts::Vector3 seg3Local = invTransform.linear() * dir;
  // Scale unit direction to the actual segemnt in the (depletion) zone
  seg3Local *= depletion / std::cos(Acts::VectorHelpers::theta(seg3Local));
  // Calulate local entry/exit before drift
  Acts::Vector2 entry = pos2Local - 0.5 * seg3Local.segment<2>(0);
  Acts::Vector2 exit = pos2Local + 0.5 * seg3Local.segment<2>(0);
  // Actually apply a drift
  // - dirftDir is assumed in local coordinates
  if (not driftDir.isApprox(Acts::Vector3(0., 0., 0.)) and
      not driftDir.isApprox(Acts::Vector3(0., 0., 1.)) and
      not driftDir.isApprox(Acts::Vector3(0., 0., -1.))) {
    // Apply the scaled drift
    auto applyDrift = [&](Acts::Vector2& local) {
      auto scaledDriftDir =
          driftDir * depletion / std::cos(Acts::VectorHelpers::theta(driftDir));
      local += scaledDriftDir.segment<2>(0);
    };

    if (driftDir.z() > 0.) {
      applyDrift(entry);
    }
    if (driftDir.z() < 0.) {
      applyDrift(exit);
    }
  }

  return {entry, exit};
}
