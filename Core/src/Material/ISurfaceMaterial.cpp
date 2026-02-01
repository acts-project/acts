// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/ISurfaceMaterial.hpp"

#include "Acts/Definitions/Common.hpp"

#include <stdexcept>

namespace Acts {

double ISurfaceMaterial::factor(Direction pDir, MaterialUpdateMode mode) const {
  if (mode == MaterialUpdateMode::NoUpdate) {
    return 0.;
  } else if (mode == MaterialUpdateMode::FullUpdate) {
    return 1.;
  } else if (mode == MaterialUpdateMode::PreUpdate) {
    return pDir == Direction::Negative() ? m_splitFactor : 1 - m_splitFactor;
  } else if (mode == MaterialUpdateMode::PostUpdate) {
    return pDir == Direction::Positive() ? m_splitFactor : 1 - m_splitFactor;
  }

  throw std::logic_error(
      "ISurfaceMaterial::factor: Unknown MaterialUpdateMode");
}

MaterialSlab ISurfaceMaterial::materialSlab(const Vector2& lp, Direction pDir,
                                            MaterialUpdateMode mode) const {
  // The plain material properties associated to this bin
  MaterialSlab plainMatProp = materialSlab(lp);
  // Scale if you have material to scale
  if (!plainMatProp.isVacuum()) {
    double scaleFactor = factor(pDir, mode);
    if (scaleFactor == 0.) {
      return MaterialSlab::Nothing();
    }
    plainMatProp.scaleThickness(scaleFactor);
  }
  return plainMatProp;
}

MaterialSlab ISurfaceMaterial::materialSlab(const Vector3& gp, Direction pDir,
                                            MaterialUpdateMode mode) const {
  // The plain material properties associated to this bin
  MaterialSlab plainMatProp = materialSlab(gp);
  // Scale if you have material to scale
  if (!plainMatProp.isVacuum()) {
    double scaleFactor = factor(pDir, mode);
    if (scaleFactor == 0.) {
      return MaterialSlab::Nothing();
    }
    plainMatProp.scaleThickness(scaleFactor);
  }
  return plainMatProp;
}

}  // namespace Acts
