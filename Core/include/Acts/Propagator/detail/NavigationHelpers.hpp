// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <variant>

namespace Acts::detail {

using AnyIntersectionObject =
    std::variant<const Surface*, const Layer*, const BoundarySurface*>;

/// @brief A candidate object for navigation
struct NavigationObjectCandidate {
  AnyIntersectionObject object;
  const Surface* representation = nullptr;
  BoundaryTolerance boundaryTolerance;

  NavigationObjectCandidate(const Surface& _surface,
                            BoundaryTolerance _boundaryTolerance)
      : object(&_surface),
        representation(&_surface),
        boundaryTolerance(_boundaryTolerance) {}
  NavigationObjectCandidate(const Layer& _layer, const Surface& _representation,
                            BoundaryTolerance _boundaryTolerance)
      : object(&_layer),
        representation(&_representation),
        boundaryTolerance(_boundaryTolerance) {}
  NavigationObjectCandidate(const BoundarySurface& _boundary,
                            BoundaryTolerance _boundaryTolerance)
      : object(&_boundary),
        representation(&_boundary.surfaceRepresentation()),
        boundaryTolerance(_boundaryTolerance) {}

  MultiIntersection3D intersect(const GeometryContext& gctx,
                                const Vector3& position,
                                const Vector3& direction,
                                double tolerance) const {
    return representation->intersect(gctx, position, direction,
                                     boundaryTolerance, tolerance);
  }

  NavigationTarget target(const Intersection3D& intersection,
                          IntersectionIndex intersectionIndex) const {
    return std::visit(
        overloaded{[&](const Surface* surface) -> NavigationTarget {
                     return {intersection, intersectionIndex, *surface,
                             boundaryTolerance};
                   },
                   [&](const Layer* layer) -> NavigationTarget {
                     return {intersection, intersectionIndex, *layer,
                             *representation, boundaryTolerance};
                   },
                   [&](const BoundarySurface* boundary) -> NavigationTarget {
                     return {intersection, intersectionIndex, *boundary,
                             boundaryTolerance};
                   }},
        object);
  }
};

/// @brief Emplace all navigation candidates for a given volume
inline void emplaceAllVolumeCandidates(
    std::vector<NavigationObjectCandidate>& candidates,
    const TrackingVolume& volume, bool resolveSensitive, bool resolveMaterial,
    bool resolvePassive,
    const BoundaryTolerance& boundaryToleranceSurfaceApproach,
    const Logger& logger) {
  // Get all boundary candidates
  {
    ACTS_VERBOSE("Searching for boundaries.");

    const auto& boundaries = volume.boundarySurfaces();

    ACTS_VERBOSE("Found " << boundaries.size() << " boundaries.");

    for (const auto& boundary : boundaries) {
      candidates.emplace_back(*boundary, BoundaryTolerance::None());
    }
  }

  // Get all layer candidates
  {
    ACTS_VERBOSE("Searching for layers.");

    const auto& layers = volume.confinedLayers()->arrayObjects();

    ACTS_VERBOSE("Found " << layers.size() << " layers.");

    for (const auto& layer : layers) {
      if (!layer->resolve(resolveSensitive, resolveMaterial, resolvePassive)) {
        continue;
      }

      if (!resolveSensitive ||
          layer->surfaceRepresentation().surfaceMaterial() != nullptr) {
        candidates.emplace_back(*layer, layer->surfaceRepresentation(),
                                boundaryToleranceSurfaceApproach);
      }

      if (layer->approachDescriptor() != nullptr) {
        const auto& approaches =
            layer->approachDescriptor()->containedSurfaces();

        for (const Surface* approach : approaches) {
          candidates.emplace_back(*layer, *approach,
                                  boundaryToleranceSurfaceApproach);
        }
      }

      // Get all surface candidates from layers
      {
        ACTS_VERBOSE("Searching for surfaces.");

        if (layer->surfaceArray() != nullptr) {
          const auto& surfaces = layer->surfaceArray()->surfaces();

          ACTS_VERBOSE("Found " << surfaces.size() << " surfaces.");

          for (const Surface* surface : surfaces) {
            candidates.emplace_back(*surface, boundaryToleranceSurfaceApproach);
          }
        }
      }
    }
  }
}

}  // namespace Acts::detail
