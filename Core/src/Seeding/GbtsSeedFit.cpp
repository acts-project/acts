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
    const Vector2& fieldDirection, const MagneticFieldProvider& magneticField,
    MagneticFieldProvider::Cache& fieldCache) {
  if (!fit.valid) {
    return std::nullopt;
  }

  const Result<Vector3> field =
      magneticField.getField(fieldPosition, fieldCache);
  if (!field.ok()) {
    return std::nullopt;
  }
  const double transverseField = field->head<2>().dot(fieldDirection);
  const double bendingField = field->z() - fit.cotTheta * transverseField;
  if (bendingField == 0) {
    return std::nullopt;
  }

  const Vector3 direction =
      Vector3(std::cos(fit.phi), std::sin(fit.phi), fit.cotTheta).normalized();
  const double qOverPt = fit.curvature / bendingField;

  FreeVector params = FreeVector::Zero();
  params.segment<3>(eFreePos0) = fit.position;
  params.segment<3>(eFreeDir0) = direction;
  params[eFreeQOverP] = qOverPt / fastHypot(1., fit.cotTheta);
  return params;
}

}  // namespace Acts::Experimental
