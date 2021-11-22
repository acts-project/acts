// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"

#include <memory>
#include <vector>

namespace Acts {

class DetectorVolume;
class VolumeBlueprint;

class BlueprintDetectorBuilder {

  /// Build a list of volumes from a volume blue print
  ///
  /// @param gctx the geometry context
  /// @param vbp the volume blueprint
  ///
  /// @return the vector of detector volumes
  static std::vector<std::shared_ptr<DetectorVolume>> build(
      const GeometryContext& gctx, const VolumeBlueprint& vbp);

};

}  // namespace Acts
