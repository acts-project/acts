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
    const GbtsSeedFit& fit, const Vector3& fieldPosition,
    const MagneticFieldProvider& magneticField,
    MagneticFieldProvider::Cache& fieldCache) {
  if (!fit.valid) {
    return std::nullopt;
  }

  const Result<Vector3> field =
      magneticField.getField(fieldPosition, fieldCache);
  if (!field.ok()) {
    return std::nullopt;
  }
  const double r = fastHypot(fieldPosition.x(), fieldPosition.y());
  const double br =
      r > 0 ? field->head<2>().dot(fieldPosition.head<2>()) / r : 0.;
  const double bendingField = field->z() - fit.cotTheta * br;
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
