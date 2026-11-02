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

}  // namespace Acts
