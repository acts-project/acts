// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

namespace TGeoPrimitivesHelpers {

inline Transform3D makeTransform(const Eigen::Vector3d& rotationMatrixCol0,
                                 const Eigen::Vector3d& rotationMatrixCol1,
                                 const Eigen::Vector3d& rotationMatrixCol2,
                                 const Eigen::Vector3d& translation) {
  Transform3D trf;
  trf.matrix().block(0, 0, 3, 1) = rotationMatrixCol0;
  trf.matrix().block(0, 1, 3, 1) = rotationMatrixCol1;
  trf.matrix().block(0, 2, 3, 1) = rotationMatrixCol2;
  trf.matrix().block(0, 3, 3, 1) = translation;
  return trf;
}

}  // namespace TGeoPrimitivesHelpers
}  // namespace Acts
