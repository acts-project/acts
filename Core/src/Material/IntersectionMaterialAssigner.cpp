// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/IntersectionMaterialAssigner.hpp"

#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <algorithm>

namespace Acts {

namespace {

struct SurfaceIntersection {
  Intersection3D intersection;
  const Surface* surface;

  constexpr static bool pathLengthOrder(
      const SurfaceIntersection& aIntersection,
      const SurfaceIntersection& bIntersection) noexcept {
    return Intersection3D::pathLengthOrder(aIntersection.intersection,
                                           bIntersection.intersection);
  }
};

std::vector<SurfaceIntersection> forwardOrderedIntersections(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const std::vector<const Surface*>& surfaces) {
  // First deal with the surface intersections
  std::vector<SurfaceIntersection> surfaceIntersections;
  // Intersect the surfaces
  for (const Surface* surface : surfaces) {
    // Get the intersection
    MultiIntersection3D multiIntersection = surface->intersect(
        gctx, position, direction, BoundaryTolerance::None());

    // Take the closest
    const Intersection3D& intersection = multiIntersection.closestForward();
    if (intersection.status() >= IntersectionStatus::reachable &&
        intersection.pathLength() > 0) {
      surfaceIntersections.emplace_back(intersection, surface);
      continue;
    }
  }
  // Sort the intersection along the pathlength
  std::ranges::sort(surfaceIntersections, SurfaceIntersection::pathLengthOrder);
  return surfaceIntersections;
}

}  // namespace

std::pair<std::vector<IAssignmentFinder::SurfaceAssignment>,
          std::vector<IAssignmentFinder::VolumeAssignment>>
IntersectionMaterialAssigner::assignmentCandidates(
    const GeometryContext& gctx, const MagneticFieldContext& /*mctx*/,
    const Vector3& position, const Vector3& direction) const {
  // The resulting candidates
  std::pair<std::vector<IAssignmentFinder::SurfaceAssignment>,
            std::vector<IAssignmentFinder::VolumeAssignment>>
      candidates;

  ACTS_DEBUG("Finding material assignment from position "
             << toString(position) << " and direction " << toString(direction));

  // Try the surfaces first
  auto sIntersections =
      forwardOrderedIntersections(gctx, position, direction, m_cfg.surfaces);
  candidates.first.reserve(sIntersections.size());
  for (auto& sIntersection : sIntersections) {
    candidates.first.push_back(IAssignmentFinder::SurfaceAssignment{
        sIntersection.surface, sIntersection.intersection.position(),
        direction});
  }

  // Now deal with the volume intersections : tracking volume first
  if (!m_cfg.trackingVolumes.empty()) {
    for (auto& trackingVolume : m_cfg.trackingVolumes) {
      // Collect the boundary surfaces
      auto boundarySurfaces = trackingVolume->boundarySurfaces();
      std::vector<const Surface*> tSurfaces;
      for (auto& boundarySurface : boundarySurfaces) {
        tSurfaces.push_back(&(boundarySurface->surfaceRepresentation()));
      }
      // Get the intersections
      auto tIntersections =
          forwardOrderedIntersections(gctx, position, direction, tSurfaces);
      // Entry/exit exists in forward direction
      if (tIntersections.size() == 2u) {
        candidates.second.push_back(IAssignmentFinder::VolumeAssignment{
            InteractionVolume(trackingVolume),
            tIntersections[0u].intersection.position(),
            tIntersections[1u].intersection.position()});
      }
    }
  }

  ACTS_DEBUG("Found " << candidates.first.size() << " surface candidates and "
                      << candidates.second.size() << " volume candidates");

  // Return the result
  return candidates;
}

}  // namespace Acts
