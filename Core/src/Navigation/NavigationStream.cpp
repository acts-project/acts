// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/NavigationStream.hpp"

#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <algorithm>
#include <unordered_set>

namespace Acts {

bool NavigationStream::initialize(const GeometryContext& gctx,
                                  const QueryPoint& queryPoint,
                                  const BoundaryTolerance& cTolerance,
                                  const double onSurfaceTolerance) {
  // Position and direction from the query point
  const Vector3& position = queryPoint.position;
  const Vector3& direction = queryPoint.direction;

  // A container collecting additional candidates from multiple
  // valid intersections
  std::vector<NavigationTarget> additionalCandidates = {};
  additionalCandidates.reserve(m_candidates.size());
  std::unordered_set<const Surface*> processed{};
  for (auto& candidate : m_candidates) {
    // Get the surface from the object intersection
    const Surface& surface = candidate.surface();
    // Check whether the surface already has been processed
    if (!processed.insert(&surface).second) {
      continue;
    }
    // Intersect the surface
    auto multiIntersection = surface.intersect(gctx, position, direction,
                                               cTolerance, onSurfaceTolerance);

    bool firstValid = multiIntersection.at(0).isValid();
    bool secondValid = multiIntersection.at(1).isValid();
    if (firstValid && !secondValid) {
      if (multiIntersection.at(0).pathLength() < -onSurfaceTolerance) {
        continue;
      }
      candidate.intersection() = multiIntersection.at(0);
      candidate.intersectionIndex() = 0;
    } else if (!firstValid && secondValid) {
      if (multiIntersection.at(1).pathLength() < -onSurfaceTolerance) {
        continue;
      }
      candidate.intersection() = multiIntersection.at(1);
      candidate.intersectionIndex() = 1;
    } else {
      // Split them into valid intersections, keep track of potentially
      // additional candidates
      bool originalCandidateUpdated = false;
      for (auto [intersectionIndex, intersection] :
           enumerate(multiIntersection)) {
        // Skip negative solutions, respecting the on surface tolerance
        if (intersection.pathLength() < -onSurfaceTolerance) {
          continue;
        }
        // Valid solution is either on surface or updates the distance
        if (intersection.isValid()) {
          if (!originalCandidateUpdated) {
            candidate.intersection() = intersection;
            candidate.intersectionIndex() = intersectionIndex;
            originalCandidateUpdated = true;
          } else {
            NavigationTarget additionalCandidate = candidate;
            additionalCandidate.intersection() = intersection;
            additionalCandidate.intersectionIndex() = intersectionIndex;
            additionalCandidates.emplace_back(additionalCandidate);
          }
        }
      }
    }
  }

  // Append the multi intersection candidates
  m_candidates.insert(m_candidates.end(), additionalCandidates.begin(),
                      additionalCandidates.end());

  // Sort the candidates by path length
  std::ranges::sort(m_candidates, NavigationTarget::pathLengthOrder);

  // If we have duplicates, we expect them to be close by in path length, so we
  // don't need to re-sort Remove duplicates on basis of the surface pointer

  /// But but but... What about the surfaces with multiple intersections?
  auto nonUniqueRange = std::ranges::unique(
      m_candidates.begin(), m_candidates.end(),
      [](const NavigationTarget& a, const NavigationTarget& b) {
        return &a.surface() == &b.surface();
      });
  m_candidates.erase(nonUniqueRange.begin(), nonUniqueRange.end());

  // The we find the first invalid candidate
  auto firstInvalid = std::ranges::find_if(
      m_candidates,
      [](const NavigationTarget& a) { return !a.intersection().isValid(); });

  // Set the range and initialize
  m_candidates.resize(std::distance(m_candidates.begin(), firstInvalid),
                      NavigationTarget::None());

  m_currentIndex = 0;
  if (m_candidates.empty()) {
    return false;
  }
  return true;
}

bool NavigationStream::update(const GeometryContext& gctx,
                              const QueryPoint& queryPoint,
                              double onSurfaceTolerance) {
  // Loop over the (currently valid) candidates and update
  for (; m_currentIndex < m_candidates.size(); ++m_currentIndex) {
    // Get the candidate, and resolve the tuple
    NavigationTarget& candidate = currentCandidate();
    // Get the surface from the object intersection
    const Surface& surface = candidate.surface();
    // (re-)Intersect the surface
    auto multiIntersection =
        surface.intersect(gctx, queryPoint.position, queryPoint.direction,
                          candidate.boundaryTolerance(), onSurfaceTolerance);
    // Split them into valid intersections
    for (auto [intersectionIndex, intersection] :
         enumerate(multiIntersection)) {
      // Skip wrong index solution
      if (intersectionIndex != candidate.intersectionIndex()) {
        continue;
      }
      // Valid solution is either on surface or updates the distance
      if (intersection.isValid()) {
        candidate.intersection() = intersection;
        return true;
      }
    }
  }
  // No candidate was reachable
  return false;
}

void NavigationStream::reset() {
  m_candidates.clear();
  m_currentIndex = 0;
}

void NavigationStream::addSurfaceCandidate(
    const Surface& surface, const BoundaryTolerance& bTolerance) {
  m_candidates.emplace_back(Intersection3D::Invalid(), 0, surface, bTolerance);
}

void NavigationStream::addSurfaceCandidates(
    std::span<const Surface*> surfaces, const BoundaryTolerance& bTolerance) {
  m_candidates.reserve(m_candidates.size() + surfaces.size());
  std::ranges::for_each(surfaces, [&](const Surface* surface) {
    m_candidates.emplace_back(Intersection3D::Invalid(), 0, *surface,
                              bTolerance);
  });
}

void NavigationStream::addPortalCandidate(const Portal& portal) {
  m_candidates.emplace_back(Intersection3D::Invalid(), 0, portal,
                            BoundaryTolerance::None());
}

AppendOnlyNavigationStream::AppendOnlyNavigationStream(NavigationStream& stream)
    : m_stream{&stream} {}

void AppendOnlyNavigationStream::addPortalCandidate(const Portal& portal) {
  m_stream->addPortalCandidate(portal);
}

void AppendOnlyNavigationStream::addSurfaceCandidate(
    const Surface& surface, const BoundaryTolerance& bTolerance) {
  m_stream->addSurfaceCandidate(surface, bTolerance);
}

}  // namespace Acts
