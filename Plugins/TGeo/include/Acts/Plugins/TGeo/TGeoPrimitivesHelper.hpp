// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <vector>
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

struct TGeoPrimitivesHelper {
  /// Helper method to create a transform from
  /// Rotation matrix vectors:
  /// @param rotationMatrixCol0
  /// @param rotationMatrixCol1
  /// @param rotationMatrixCol2
  /// And translation
  /// @param translation
  static inline Transform3D makeTransform(
      const Eigen::Vector3d& rotationMatrixCol0,
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

  /// Private helper method : match string with wildcards
  /// @param first is the one with the potential wildcard
  /// @param second is the test string
  static bool match(const char* first, const char* second);

  /// Private helper method : match string with wildcards
  /// Method that uses the match method with wild cards and
  /// performs it on an input list
  /// @param first is the one with the potential wildcard
  /// @param second is the test string
  static bool match(const std::vector<std::string>& first, const char* second);
};
}  // namespace Acts
