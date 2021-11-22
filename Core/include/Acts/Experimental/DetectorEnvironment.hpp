// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Experimental/DetectorEnvironment.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <vector>
#include <limits>

namespace Acts {

class DetectorVolume;
class Portal;
class Surface;

using SurfaceIntersection = ObjectIntersection<Surface>;
using PortalIntersection = ObjectIntersection<Portal, Surface>;

/// A pure navigation struct, that describes the current
/// environment, it provided by the portal at entry into into
/// detector volume.
struct DetectorEnvironment {
  
  /// The status of the environment
  enum Status : unsigned int {
    eUninitialized = 0,
    eTowardsSurface = 1,
    eOnSurface = 2,
    eTowardsPortal = 3,
    eOnPortal = 4
  };

  /// That is the volume you are processing
  const DetectorVolume* volume = nullptr;

  /// That are the candidate surfaces to process
  std::vector<SurfaceIntersection> surfaces = {};

  /// That are the portals for leaving that environment
  std::vector<PortalIntersection> portals = {};

  /// Set the current surface
  const Surface* currentSurface = nullptr;

  /// Indicate the status of this environment
  Status status = eUninitialized;

};

}  // namespace Acts
