// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"

#include <any>
#include <cstddef>
#include <vector>

namespace Acts {

class Surface;

namespace Experimental {

class Portal;
class Detector;
class DetectorVolume;

/// @brief A navigation state struct that is holding the current navigation information
///
/// It relies on Surfaces and Portals, all navigation entities have to be
/// described in these terms.
struct NavigationState {
  /// Surface candidate vector alias, this allows to use e.g. boost_small vector
  /// or other stl like containers
  using SurfaceCandidates = std::vector<NavigationTarget>;

  /// The current position
  Vector3 position = Vector3(0., 0., 0.);

  /// The current direction
  Vector3 direction = Vector3(0., 0., 0.);

  /// The current detector in processing
  const Detector* currentDetector = nullptr;

  /// The current volume in processing, i.e. the position is inside
  const DetectorVolume* currentVolume = nullptr;

  /// The current surface, i.e the position is on surface
  const Surface* currentSurface = nullptr;

  /// The current portal, i.e the position is on portal
  const Portal* currentPortal = nullptr;

  /// That are the candidate surfaces to process
  SurfaceCandidates surfaceCandidates = {};

  /// Starting index of the surface candidate to correctly identify the first
  /// surface
  int surfaceCandidateIndex = -1;

  /// Boundary directives for surfaces
  BoundaryTolerance surfaceBoundaryTolerance = BoundaryTolerance::None();

  /// An overstep tolerance
  double overstepTolerance = -100 * UnitConstants::um;

  /// Auxiliary attached information
  std::any auxiliary;

  /// Get the current surface candidate being processed
  /// @return Reference to the current surface candidate
  const NavigationTarget& surfaceCandidate() const {
    return surfaceCandidates.at(surfaceCandidateIndex);
  }
};

}  // namespace Experimental
}  // namespace Acts
