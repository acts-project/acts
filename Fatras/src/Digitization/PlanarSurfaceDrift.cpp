// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Digitization/PlanarSurfaceDrift.hpp"

#include "Acts/Utilities/Helpers.hpp"
#include "ActsFatras/Digitization/DigitizationError.hpp"

#include <cmath>

Acts::Result<std::tuple<ActsFatras::PlanarSurfaceDrift::Segment2D,
                        ActsFatras::PlanarSurfaceDrift::Segment3D>>
ActsFatras::PlanarSurfaceDrift::toReadout(const Acts::GeometryContext& gctx,
                                          const Acts::Surface& surface,
                                          double thickness,
                                          const Acts::Vector3& pos,
                                          const Acts::Vector3& dir,
                                          const Acts::Vector3& driftDir) const {
  // Transform the hit & direction into the local surface frame
  const auto& invTransform = surface.transform(gctx).inverse();
  Acts::Vector3 pos3Local(invTransform * pos);
  Acts::Vector3 seg3Local = invTransform.linear() * dir.normalized();
  if (std::abs(seg3Local.z()) < Acts::s_epsilon) {
    // The segment is parallel to the surface
    return ActsFatras::DigitizationError::DriftError;
  }
  // Scale the unit vector to the thickness of the module
  double scale = thickness / seg3Local.z();
  seg3Local *= scale;
  // The drift direction is in the local frame, so we need to transform it
  Acts::Vector3 entry = pos3Local - 0.5 * seg3Local;
  Acts::Vector3 exit = pos3Local + 0.5 * seg3Local;
  Acts::Vector3 driftedEntry = entry;
  Acts::Vector3 driftedExit = exit;
  // Apply the drift if it is not zero
  if (Acts::VectorHelpers::perp(driftDir) > Acts::s_epsilon &&
      std::abs(driftDir.z()) > Acts::s_epsilon) {
    // Apply the drift to the entry and exit points
    double driftScale = 0.5 * thickness / driftDir.z();
    driftedEntry += driftScale * driftDir;
    driftedExit -= driftScale * driftDir;
  }
  // The drifted segment and the original segment
  return std::tuple{
      Segment2D{driftedEntry.segment<2>(0), driftedExit.segment<2>(0)},
      Segment3D{entry, exit}};
}
