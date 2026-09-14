// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/MagneticFieldProvider.hpp"

namespace Acts {

Result<MagneticFieldProvider::FieldAndGradient>
MagneticFieldProvider::getFieldAndGradient(const Vector3& /*position*/,
                                           Cache& /*cache*/) const {
  return Result<FieldAndGradient>::failure(MagneticFieldError::NotImplemented);
}

Result<MagneticFieldProvider::FieldAndGradient> getFieldAndGradientNumerically(
    const MagneticFieldProvider& provider, const Vector3& position,
    MagneticFieldProvider::Cache& cache, double epsilon) {
  using FieldAndGradient = MagneticFieldProvider::FieldAndGradient;

  FieldAndGradient result;

  Result<Vector3> field = provider.getField(position, cache);
  if (!field.ok()) {
    return Result<FieldAndGradient>::failure(field.error());
  }
  result.field = *field;

  for (std::size_t j = 0; j < 3; ++j) {
    Vector3 delta = Vector3::Zero();
    delta[j] = epsilon;

    Result<Vector3> plus = provider.getField(position + delta, cache);
    if (!plus.ok()) {
      return Result<FieldAndGradient>::failure(plus.error());
    }
    Result<Vector3> minus = provider.getField(position - delta, cache);
    if (!minus.ok()) {
      return Result<FieldAndGradient>::failure(minus.error());
    }

    result.gradient.col(j) = (*plus - *minus) / (2. * epsilon);
  }

  return Result<FieldAndGradient>::success(result);
}

}  // namespace Acts
