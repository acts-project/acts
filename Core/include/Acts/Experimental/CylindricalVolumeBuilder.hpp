// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include <memory>
#include <string>

/// @note This file is foreseen for the `Geometry` module

namespace Acts {

class DetectorVolume;
class LayerBlueprint;

struct CylindricalVolumeBuilder {

  /// Helper function that creates a cylindrical volumes from ProtoLayers
  /// 
  /// @param protoLayers the protoLayers to be transformed into volumes
  /// @param name the base name of the volumes
  ///
  /// @return the new cylindrical volumes
  static std::vector<std::shared_ptr<DetectorVolume>> volumesInZ(
      const std::vector<LayerBlueprint>& protoLayers,      
      const std::string& name = "Volume");

  /// Helper function that creates a cylindrical container ordered in R
  /// 
  /// @param protoLayers the protoLayers to be transformed into volumes
  /// @param name the base name of the volumes
  ///
  /// @return the new (cylindrical) container volume
  static std::vector<std::shared_ptr<DetectorVolume>> volumesInR(
      const std::vector<LayerBlueprint>& protoLayers,      
      const std::string& name = "Volume");

};

}  // namespace Acts