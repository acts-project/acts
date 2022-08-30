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

/// Declare a navigation state updator
///
/// This delegate dispatches the local navigation action
/// to a dedicated struct or function that is optimised for
/// the given environment.
///
/// @param nState is the navigation state to be updated
/// @param volume is the volume for which this should be called
/// @param gctx is the current geometry context
/// @param position is the position at the query
/// @param direction is the direction at the query
/// @param absMomentum is the absolute momentum at query
/// @param charge is the charge to be used for the intersection
///
using NavigationStateUpdator = Delegate<void(
    NavigationState& nState, const DetectorVolume& volume,
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, ActsScalar absMomentum, ActsScalar charge)>;

using NavigationStateUpdatorStore = std::shared_ptr<void>;

/// Declare a Detctor Volume Switching delegate
///
/// @param gctx is the current geometry context
/// @param position is the position at the query
/// @param direction is the direction at the query
///
/// @return the new DetectorVolume into which one changes at this switch
using DetectorVolumeLink = Delegate<const DetectorVolume*(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction)>;

/// @brief Definition of the link store for ownership control
using DetectorVolumeLinkStore = std::shared_ptr<void>;

}  // namespace Experimental
}  // namespace Acts
