// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {
namespace Experimental {

/// A current detector building snapshot
///
/// It captures the current built volumes and their
/// currently set outside detector skin
///
struct BlueprintBlock {
  std::vector<std::shared_ptr<DetectorVolume>> volumes = {};
  std::vector<std::shared_ptr<Portal>> skin = {};

  /// Build a simple shell from a volume
  /// @param dVolume the detector volume
  BlueprintBlock(DetectorVolume& dVolume)
      : volumes({dVolume.getSharedPtr()}), skin(dVolume.portalPtrs()) {}

  // Empty constructor
  BlueprintBlock() = default;
};

}  // namespace Experimental
}  // namespace Acts
