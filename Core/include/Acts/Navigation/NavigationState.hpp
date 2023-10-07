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
  /// @brief  A surface candidate and its intersection
  ///
  /// A candidates can either be a surface or a portal (which contain a surface)
  struct SurfaceCandidate {
    /// A candidate intersection, in Surface view
    ObjectIntersection<Surface> objectIntersection;
    /// A candidate is either a detector Surface
    const Surface* surface = nullptr;
    /// Or a portal
    const Portal* portal = nullptr;
    /// The boundary check used for the candidate, boundary checks
    /// can differ for sensitive surfaces and portals
    BoundaryCheck boundaryCheck = true;
  };

  /// Surface candidate vector alias, this allows to use e.g. boost_small vector
  /// or other stl like containers
  using SurfaceCandidates = std::vector<SurfaceCandidate>;

  /// The current position
  Vector3 position = Vector3(0., 0., 0.);

  /// The current direction
  Vector3 direction = Vector3(0., 0., 0.);

  /// The current absolute momentum
  ActsScalar absMomentum = 0.;

  /// The current charge
  ActsScalar charge = 0.;

  /// The current magnetic field
  Vector3 magneticField = Vector3(0., 0., 0.);

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
