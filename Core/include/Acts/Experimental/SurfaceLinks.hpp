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
  /// @param maxPathLength the maximal path length
  /// @note provideAll is (obvisouly) ignored in this context
  ///
  /// @note the onSurface solution will be given first
  ///
  /// @return the surface intersections (ordered)
  std::vector<SurfaceIntersection> operator()(
      const GeometryContext& gctx, const DetectorVolume& volume,
      const Vector3& position, const Vector3& direction,
      const BoundaryCheck& bCheck, const ActsScalar maxPathLength, bool) const {
    // A volume needs to be present and so do surfaces
    if (not volume.surfaces().empty()) {
      std::vector<SurfaceIntersection> sIntersections;
      sIntersections.reserve(guessedNumberOfCandidates);
      // Run over all surfaces and create the trial & error candidate lists in
      for (const auto& s : volume.surfaces()) {
        auto sIntersection = s->intersect(gctx, position, direction, bCheck);
        // The intersection needs to be valid and within range
        if (sIntersection.intersection.status >=
                Intersection3D::Status::reachable and
            sIntersection.intersection.pathLength < maxPathLength) {
          sIntersections.push_back(sIntersection);
        }
      }
      // Sort and return
      std::sort(sIntersections.begin(), sIntersections.end());
      return sIntersections;
    }
    return {};
  }
};

}  // namespace Acts
