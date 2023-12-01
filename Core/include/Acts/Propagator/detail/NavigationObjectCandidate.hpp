// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/AnyIntersection.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"

#include <variant>

namespace Acts {
namespace detail {

/// @brief A candidate object for navigation
///
/// TODO if the intersection would carry a `std::any` we would not need this
/// class
struct NavigationObjectCandidate {
  using SurfaceObject = const Surface*;
  using LayerObject = const Layer*;
  using BoundaryObject = const BoundarySurfaceT<TrackingVolume>*;
  using AnyObject = std::variant<SurfaceObject, LayerObject, BoundaryObject>;

  AnyObject object;
  const Surface* representation = nullptr;
  BoundaryCheck boundaryCheck;

  NavigationObjectCandidate(AnyObject _object, const Surface* _representation,
                            BoundaryCheck _boundaryCheck)
      : object(_object),
        representation(_representation),
        boundaryCheck(std::move(_boundaryCheck)) {}

  AnyMultiIntersection intersect(const GeometryContext& gctx,
                                 const Vector3& position,
                                 const Vector3& direction,
                                 ActsScalar tolerance) const {
    if (std::holds_alternative<SurfaceObject>(object)) {
      const auto& surface = std::get<SurfaceObject>(object);
      auto intersection = representation->intersect(gctx, position, direction,
                                                    boundaryCheck, tolerance);
      return AnyMultiIntersection(SurfaceMultiIntersection(
          intersection.intersections(), surface, representation));
    }
    if (std::holds_alternative<LayerObject>(object)) {
      const auto& layer = std::get<LayerObject>(object);
      auto intersection = representation->intersect(gctx, position, direction,
                                                    boundaryCheck, tolerance);
      return AnyMultiIntersection(LayerMultiIntersection(
          intersection.intersections(), layer, representation));
    }
    if (std::holds_alternative<BoundaryObject>(object)) {
      const auto& boundary = std::get<BoundaryObject>(object);
      auto intersection = representation->intersect(gctx, position, direction,
                                                    boundaryCheck, tolerance);
      return AnyMultiIntersection(BoundaryMultiIntersection(
          intersection.intersections(), boundary, representation));
    }
    throw std::runtime_error("unknown type");
  }
};

/// @brief Emplace all navigation candidates for a given volume
inline void emplaceAllVolumeCandidates(
    std::vector<detail::NavigationObjectCandidate>& candidates,
    const TrackingVolume& volume, bool resolveSensitive, bool resolveMaterial,
    bool resolvePassive, BoundaryCheck boundaryCheckSurfaceApproach,
    const Logger& logger) {
  auto addCandidate = [&](detail::NavigationObjectCandidate::AnyObject object,
                          const Surface* representation,
                          BoundaryCheck boundaryCheck) {
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

}  // namespace detail
}  // namespace Acts
