// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/PlanarSurfaceDrift.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <cmath>

std::tuple<ActsFatras::PlanarSurfaceDrift::Segment2D, double>
ActsFatras::PlanarSurfaceDrift::toReadout(const Acts::GeometryContext& gctx,
                                          const Acts::Surface& surface,
                                          double thickness,
                                          const Acts::Vector3& pos,
                                          const Acts::Vector3& dir,
                                          [[ maybe_unused ]] const Acts::Vector3& driftDir) const {
  // Transform the hit & direction into the local surface frame
  const auto& invTransform = surface.transform(gctx).inverse();
  Acts::Vector2 pos2Local = (invTransform * pos).segment<2>(0);
  Acts::Vector3 seg3Local = invTransform.linear() * dir;
  double scale = thickness / seg3Local.z();
  seg3Local *= scale;

  Acts::Vector2 seg2Local = seg3Local.segment<2>(0);
  Acts::Vector2 entry = pos2Local - 0.5 * seg3Local.segment<2>(0);
  Acts::Vector2 exit = pos2Local + 0.5 * seg3Local.segment<2>(0);
  double driftScale = seg3Local.norm() / seg2Local.norm();

  return { ActsFatras::PlanarSurfaceDrift::Segment2D{entry, exit}, driftScale };
}
