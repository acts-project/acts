// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Experimental/InternalBlueprint.hpp"
#include "Acts/Experimental/SurfaceLinks.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include <vector>
#include <memory>

namespace Acts {

class Surface;
class VolumeBounds;

struct BarrelGridSurfaceLinksGenerator {
  /// Single or multiple grid structure
  bool singleGrid = true;

  /// The number of bins in phi 
  size_t binsPhi = 1;

  /// The number of bins in z 
  size_t binsZ   = 1;

  /// This generates a grid based surface link ensemble for cylindrical
  /// volumes
  ///
  /// @param gctx is the geometry context during generation call
  /// @param surfaces are the internal surfaces in the order
  ///        as they'll appear in the volume.surfaces() call
  /// @param volumeBounds are the volume bounds of this volume
  InternalSurfaceLinks operator()(
      const GeometryContext& gctx,
      const std::vector<std::shared_ptr<Surface>>& surfaces,
      const VolumeBounds* volumeBounds) const noexcept(false);
};



}  // namespace Acts