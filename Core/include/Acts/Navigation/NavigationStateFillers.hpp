// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <vector>

namespace Acts {

class Surface;

namespace Experimental {

class Portal;
class Detector;
class DetectorVolume;

/// Filler of the current volume
struct DetectorVolumeFiller {
  /// Helper struct that allows to fill a volume into the
  /// navigation state, it allows to use common navigation
  /// structs for volume, portal, surfaces
  ///
  /// @param nState the navigation state
  /// @param volume the volume that is filled
  inline static void fill(NavigationState& nState,
                          const DetectorVolume* volume) {
    nState.currentVolume = volume;
  }
};

/// Fillers and attachers for surfaces to act on the navigation state
struct SurfacesFiller {
  /// Helper struct that allows to fill surfaces into the candidate vector it
  /// allows to use common navigation structs for volume, portal, surfaces
  ///
  /// @param nState the navigation state
  /// @param surfaces the surfaces that are filled in
  inline static void fill(NavigationState& nState,
                          const std::vector<const Surface*>& surfaces) {
    std::for_each(surfaces.begin(), surfaces.end(), [&](const auto& s) {
      nState.surfaceCandidates.push_back(NavigationState::SurfaceCandidate{
          ObjectIntersection<Surface>{}, s, nullptr,
          nState.surfaceBoundaryCheck});
    });
  }
};

/// Fillers and attachers for portals to act on the navigation state
struct PortalsFiller {
  /// Helper struct that allows to fill surfaces into the candidate vector it
  /// allows to use common navigation structs for volume, portal, surfaces
  ///
  /// @param nState the navigation state
  /// @param portals the portals that are filled in
  inline static void fill(NavigationState& nState,
                          const std::vector<const Portal*>& portals) {
    std::for_each(portals.begin(), portals.end(), [&](const auto& p) {
      nState.surfaceCandidates.push_back(NavigationState::SurfaceCandidate{
          ObjectIntersection<Surface>{}, nullptr, p, true});
    });
  }
};

}  // namespace Experimental
}  // namespace Acts
