// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <string>
#include <vector>

namespace ActsPlugins {
/// @addtogroup root_plugin
/// @{

/// Helper struct for creating ACTS primitives from TGeo geometry
struct TGeoPrimitivesHelper {
  /// Helper method to create a transform from
  /// Rotation matrix vectors:
  /// @param rotationMatrixCol0
  /// @param rotationMatrixCol1
  /// @param rotationMatrixCol2
  /// And translation
  /// @param translation
  /// @return The constructed 3D transformation
  static inline Acts::Transform3 makeTransform(
      const Eigen::Vector3d& rotationMatrixCol0,
      const Eigen::Vector3d& rotationMatrixCol1,
      const Eigen::Vector3d& rotationMatrixCol2,
      const Eigen::Vector3d& translation) {
    Acts::Transform3 trf;
    trf.matrix().block<3, 1>(0, 0) = rotationMatrixCol0;
    trf.matrix().block<3, 1>(0, 1) = rotationMatrixCol1;
    trf.matrix().block<3, 1>(0, 2) = rotationMatrixCol2;
    trf.matrix().block<3, 1>(0, 3) = translation;
    return trf;
  }

  /// Private helper method : match string with wildcards
  /// @param first is the one with the potential wildcard
  /// @param second is the test string
  /// @return True if strings match (accounting for wildcards)
  static bool match(const char* first, const char* second);

  /// Private helper method : match string with wildcards
  /// Method that uses the match method with wild cards and
  /// performs it on an input list
  /// @param first is the one with the potential wildcard
  /// @param second is the test string
  /// @return True if any pattern in the list matches the test string
  static bool match(const std::vector<std::string>& first, const char* second);
};
/// @}
}  // namespace ActsPlugins
