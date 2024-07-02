// This file is part of the Acts project.
//
// Copyright (C) 2023-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <variant>

namespace Acts::detail {

using AnyIntersectionObject =
    std::variant<const Surface*, const Layer*, const BoundarySurface*>;

/// @brief A candidate object for navigation
struct NavigationObjectCandidate {
  AnyIntersectionObject object;
  const Surface* representation = nullptr;
  BoundaryCheck boundaryCheck;

  NavigationObjectCandidate(AnyIntersectionObject _object,
                            const Surface* _representation,
                            BoundaryCheck _boundaryCheck)
      : object(_object),
        representation(_representation),
        boundaryCheck(std::move(_boundaryCheck)) {}

  std::pair<SurfaceMultiIntersection, AnyIntersectionObject> intersect(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction, ActsScalar tolerance) const {
    if (std::holds_alternative<const Surface*>(object)) {
      const auto& surface = std::get<const Surface*>(object);
      auto intersection = representation->intersect(gctx, position, direction,
                                                    boundaryCheck, tolerance);
      return {intersection, surface};
    }
    if (std::holds_alternative<const Layer*>(object)) {
      const auto& layer = std::get<const Layer*>(object);
      auto intersection = representation->intersect(gctx, position, direction,
                                                    boundaryCheck, tolerance);
      return {intersection, layer};
    }
    if (std::holds_alternative<const BoundarySurface*>(object)) {
      const auto& boundary = std::get<const BoundarySurface*>(object);
      auto intersection = representation->intersect(gctx, position, direction,
                                                    boundaryCheck, tolerance);
      return {intersection, boundary};
    }
    throw std::runtime_error("unknown type");
  }
};

/// Composes an intersection and a bounds check into a navigation candidate.
/// This is used to consistently update intersections after creation.
struct IntersectionCandidate {
  SurfaceIntersection intersection;
  detail::AnyIntersectionObject anyObject;
  BoundaryCheck boundaryCheck;

  IntersectionCandidate(SurfaceIntersection _intersection,
                        detail::AnyIntersectionObject _anyObject,
                        BoundaryCheck _boundaryCheck)
      : intersection(std::move(_intersection)),
        anyObject(_anyObject),
        boundaryCheck(std::move(_boundaryCheck)) {}

  template <typename object_t>
  bool checkType() const {
    return std::holds_alternative<const object_t*>(anyObject);
  }

  template <typename object_t>
  const object_t* object() const {
    return std::get<const object_t*>(anyObject);
  }

  static bool forwardOrder(const IntersectionCandidate& aCandidate,
                           const IntersectionCandidate& bCandidate) {
    return Intersection3D::pathLengthOrder(
        aCandidate.intersection.intersection(),
        bCandidate.intersection.intersection());
  }
};

/// @brief Emplace all navigation candidates for a given volume
inline void emplaceAllVolumeCandidates(
    std::vector<detail::NavigationObjectCandidate>& candidates,
    const TrackingVolume& volume, bool resolveSensitive, bool resolveMaterial,
    bool resolvePassive, const BoundaryCheck& boundaryCheckSurfaceApproach,
    const Logger& logger) {
  auto addCandidate = [&](AnyIntersectionObject object,
                          const Surface* representation,
                          const BoundaryCheck& boundaryCheck) {
    candidates.emplace_back(object, representation, boundaryCheck);
  };

  // Get all boundary candidates
  {
    ACTS_VERBOSE("Searching for boundaries.");

    const auto& boundaries = volume.boundarySurfaces();

    ACTS_VERBOSE("Found " << boundaries.size() << " boundaries.");

    for (const auto& boundary : boundaries) {
      addCandidate(boundary.get(), &boundary->surfaceRepresentation(),
                   BoundaryCheck(true));
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
        addCandidate(layer.get(), &layer->surfaceRepresentation(),
                     BoundaryCheck(true));
      }

      if (layer->approachDescriptor() != nullptr) {
        const auto& approaches =
            layer->approachDescriptor()->containedSurfaces();

        for (const auto& approach : approaches) {
          addCandidate(layer.get(), approach, BoundaryCheck(true));
        }
      }

      // Get all surface candidates from layers
      {
        ACTS_VERBOSE("Searching for surfaces.");

        if (layer->surfaceArray() != nullptr) {
          const auto& surfaces = layer->surfaceArray()->surfaces();

          ACTS_VERBOSE("Found " << surfaces.size() << " surfaces.");

          for (const auto& surface : surfaces) {
            addCandidate(surface, surface, boundaryCheckSurfaceApproach);
          }
        }
      }
    }
  }
}

}  // namespace Acts::detail
