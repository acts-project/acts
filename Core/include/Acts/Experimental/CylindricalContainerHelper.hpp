// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <string>
#include <vector>

/// @note This file is foreseen for the `Geometry` module

namespace Acts {

class DetectorVolume;

struct CylindricalContainerHelper {
  /// Helper function that creates a cylindrical container ordered in along Z
  ///
  /// @param containerVolumes the contained volumes that will be moved
  /// @param containerName the name of the new container
  ///
  /// @return the new (cylindrical) container volume
  static std::shared_ptr<DetectorVolume> containerInZ(
      std::vector<std::shared_ptr<DetectorVolume>>&& containerVolumes,
      const std::string& containerName = "Container");

  /// Helper function that creates a cylindrical container ordered in R
  ///
  /// @param containerVolumes the contained volumes that will be moved
  /// @param containerName the name of the new container
  ///
  /// @return the new (cylindrical) container volume
  static std::shared_ptr<DetectorVolume> containerInR(
      std::vector<std::shared_ptr<DetectorVolume>>&& containerVolumes,
      const std::string& containerName = "Container");
};

}  // namespace Acts