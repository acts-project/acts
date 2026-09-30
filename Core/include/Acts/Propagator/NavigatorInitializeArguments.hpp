// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"

namespace Acts {

class Surface;
class TrackingVolume;

/// Per-run inputs to the navigator initialization.
///
/// A navigation state can serve more than one run, so everything that changes
/// from run to run is passed here instead of through the navigator options.
struct NavigatorInitializeArguments {
  /// Position the navigation run starts at
  Vector3 position{};
  /// Direction the navigation run starts with
  Vector3 direction{};
  /// Direction of the propagation along the trajectory
  Direction propagationDirection = Direction::Forward();

  /// Surface the navigation run starts on, if any
  const Surface* startSurface = nullptr;
  /// Surface the navigation run targets, if any
  const Surface* targetSurface = nullptr;
  /// Hint for the start volume which short-cuts the volume search, if known
  const TrackingVolume* startVolume = nullptr;
};

}  // namespace Acts
