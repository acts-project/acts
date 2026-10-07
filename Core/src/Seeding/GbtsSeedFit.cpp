// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsSeedFit.hpp"

#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cmath>

namespace Acts::Experimental {

std::optional<FreeVector> freeParametersFromGbtsSeedFit(
    const GbtsSeedFit& fit, std::span<const Vector3> spacePoints,
    const MagneticFieldProvider& magneticField,
    MagneticFieldProvider::Cache& fieldCache) {
  if (!fit.valid || spacePoints.empty()) {
    return std::nullopt;
  }

  double bendingField = 0;
  for (const Vector3& position : spacePoints) {
    const Result<Vector3> field = magneticField.getField(position, fieldCache);
    if (!field.ok()) {
      return std::nullopt;
    }
    const double r = fastHypot(position.x(), position.y());
    const double br =
        r > 0 ? (field->x() * position.x() + field->y() * position.y()) / r
              : 0.;
    bendingField += field->z() - fit.cotTheta * br;
  }
  bendingField /= static_cast<double>(spacePoints.size());
  if (bendingField == 0) {
    return std::nullopt;
  }

  const double sinTheta = 1. / fastHypot(1., fit.cotTheta);
  FreeVector free = FreeVector::Zero();
  free.segment<3>(eFreePos0) = fit.position;
  free[eFreeDir0] = std::cos(fit.phi) * sinTheta;
  free[eFreeDir1] = std::sin(fit.phi) * sinTheta;
  free[eFreeDir2] = fit.cotTheta * sinTheta;
  free[eFreeQOverP] = fit.curvature / bendingField * sinTheta;
  return free;
}

}  // namespace Acts::Experimental
