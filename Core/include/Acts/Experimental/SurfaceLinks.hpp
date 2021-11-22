// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace Acts {

/// Helper function to get the surface candidates from a list, can and should
/// be reuesed by all surface finders as it guaraqntees consistent path,
/// overstep and alternative solutions handling
///
/// @param gctx is the current geometry conbtext
/// @param surfaces is the pre-intersected surface candidates
/// @param position is the position at the query
/// @param direction is the direction at the query
/// @param bCheck is the BoundaryCheck
/// @param pathRange is the allowed path range for the intersection
/// @param guessedNumber is the guessed number of intersections for
/// reserving
///
/// @note onSurface solutions are ranked last
inline std::vector<SurfaceIntersection> surfaceCandidates(
    const GeometryContext& gctx, std::vector<const Surface*> surfaces,
    const Vector3& position, const Vector3& direction,
    const BoundaryCheck& bCheck, const std::array<ActsScalar, 2>& pathRange,
    size_t guessedNumber) {
  // Return object and reserving
  std::vector<SurfaceIntersection> sIntersections;
  sIntersections.reserve(guessedNumber);
  for (const auto& s : surfaces) {
    auto sIntersection = s->intersect(gctx, position, direction, bCheck);
    // Swap in case of two solutions and conflict with overstep tolerance
    if (sIntersection.intersection.pathLength + s_onSurfaceTolerance < pathRange[0] and
        sIntersection.alternative.pathLength + s_onSurfaceTolerance > pathRange[0] and
        sIntersection.alternative.status >= Intersection3D::Status::reachable) {
      sIntersection.swapSolutions();
    }
    // The intersection needs to be valid and within range
    if (sIntersection.intersection.status >=
            Intersection3D::Status::reachable and
        sIntersection.intersection.pathLength + s_onSurfaceTolerance >= pathRange[0] and
        sIntersection.intersection.pathLength - s_onSurfaceTolerance < pathRange[1]) {
      sIntersections.push_back(sIntersection);
    }
  }
  // Sort and return
  std::sort(sIntersections.begin(), sIntersections.end());
  return sIntersections;
}

struct AllSurfaces {
  /// A guess for the number of candidates (for vector reserve)
  size_t guessedNumberOfCandidates = 25;

  /// Call operator for this struct to provide the list of surfaces
  ///
  /// @param gctx the current geometry context
  /// @param volume the currenct volume
  /// @param position the current position for intersection start
  /// @param direction the current momentum for intersecton start
  /// @param bCheck the boundary check
  /// @param pathRange the allowed path range
  /// @note provideAll is (obvisouly) ignored in this context
  ///
  /// @note the onSurface solution will be given first
  ///
  /// @return the surface intersections (ordered)
  std::vector<SurfaceIntersection> operator()(
      const GeometryContext& gctx, const DetectorVolume& volume,
      const Vector3& position, const Vector3& direction,
      const BoundaryCheck& bCheck, const std::array<ActsScalar, 2>& pathRange,
      bool) const {

    // A volume needs to be present and so do surfaces
    if (not volume.surfaces().empty()) {
      return surfaceCandidates(gctx, volume.surfaces(), position, direction,
                               bCheck, pathRange, guessedNumberOfCandidates);
    }
    return {};
  }
};

}  // namespace Acts
