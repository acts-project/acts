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
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <any>
#include <vector>

#include <boost/container/small_vector.hpp>

/// @note this is foreseen for the 'Geometry' module

namespace Acts {

class Surface;

namespace Experimental {

class Portal;
class DetectorVolume;

/// @brief A navigation state holding the current information
/// about volume, surfaces, and portals
struct NavigationState {
  /// @brief  A surface candidate and its intersection
  struct SurfaceCandidate {
    /// A candidate intersection, in Surface view
    ObjectIntersection<Surface> objectIntersection;
    /// A candidate is either a detector Surface
    const Surface* surface = nullptr;
    /// Or a portal
    const Portal* portal = nullptr;
    /// The boundary check used for these
    BoundaryCheck bCheck = true;
  };

  using SurfaceCandidates =
      boost::container::small_vector<SurfaceCandidate, 16u>;

  /// The current volume in processing
  const DetectorVolume* currentVolume = nullptr;

  /// The current surface, i.e the track is on surface
  const Surface* currentSurface = nullptr;

  /// That are the candidate surfaces to process
  SurfaceCandidates surfaceCandidates = {};
  SurfaceCandidates::iterator surfaceCandidate = surfaceCandidates.end();

  /// Boundary directives for surfaces
  BoundaryCheck surfaceBoundaryCheck = true;

  /// An overstep tolerance
  ActsScalar overstepTolerance = -0.1;

  /// Auxilliary attached information
  std::any auxilliary;
};

}  // namespace Experimental
}  // namespace Acts
