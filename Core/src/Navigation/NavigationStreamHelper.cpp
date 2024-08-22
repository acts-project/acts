// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/NavigationStreamHelper.hpp"

#include "Acts/Surfaces/Surface.hpp"

bool Acts::NavigationStreamHelper::initializeStream(
    NavigationStream& stream, const GeometryContext& gctx,
    const NavigationStream::QueryPoint& queryPoint,
    BoundaryTolerance cTolerance) {
  // Position and direction from the query point
  const Vector3& position = queryPoint.position;
  const Vector3& direction = queryPoint.direction;

  // De-duplicate first (necessary to deal correctly with multiple
  // intersections) - sort them by surface pointer
  std::sort(stream.candidates.begin(), stream.candidates.end(),
            [](const NavigationStream::Candidate& a,
               const NavigationStream::Candidate& b) {
              return (&a.surface()) < (&b.surface());
            });
  // Remove duplicates on basis of the surface pointer
  stream.candidates.erase(
      std::unique(stream.candidates.begin(), stream.candidates.end(),
                  [](const NavigationStream::Candidate& a,
                     const NavigationStream::Candidate& b) {
                    return (&a.surface()) == (&b.surface());
                  }),
      stream.candidates.end());

  // A container collecting additional valid ones from multi intersections
  std::vector<NavigationStream::Candidate> miCandidates = {};
  for (auto& [sIntersection, portal, bTolerance, abortTarget] :
       stream.candidates) {
    // Get the surface from the object intersection
    const Surface* surface = sIntersection.object();
    // Intersect the surface
    auto multiIntersection = surface->intersect(
        gctx, position, direction, cTolerance, s_onSurfaceTolerance);

    // Split them into valid intersections
    bool miCandidate = false;
    for (auto& rsIntersection : multiIntersection.split()) {
      // Skip negative solutions
      if (rsIntersection.pathLength() < 0.) {
        continue;
      }
      // Valid solution is either on surface or updates the distance
      if (rsIntersection.isValid()) {
        if (!miCandidate) {
          sIntersection = rsIntersection;
          miCandidate = true;
        } else {
          miCandidates.push_back(
              NavigationStream::Candidate{rsIntersection, portal, bTolerance});
        }
      }
    }
  }

  // Append the multi intersection candidates
  stream.candidates.insert(stream.candidates.end(), miCandidates.begin(),
                           miCandidates.end());

  // Sort the candidates by path length
  std::sort(stream.candidates.begin(), stream.candidates.end(),
            [](const NavigationStream::Candidate& a,
               const NavigationStream::Candidate& b) {
              return a.intersection.pathLength() < b.intersection.pathLength();
            });

  // The we find the first invalid candidate
  auto firstInvalid = std::find_if(
      stream.candidates.begin(), stream.candidates.end(),
      [](const NavigationStream::Candidate& a) {
        const auto& [aIntersection, aPortal, aTolerance, abortTarget] = a;
        return !aIntersection.isValid();
      });

  // Set the range and initialize
  stream.candidates.resize(
      std::distance(stream.candidates.begin(), firstInvalid));

  stream.currentIndex = 0;
  if (stream.candidates.empty()) {
    return false;
  }
  return true;
}

bool Acts::NavigationStreamHelper::updateStream(
    NavigationStream& stream, const GeometryContext& gctx,
    const NavigationStream::QueryPoint& queryPoint) {
  // Position and direction from the query point
  const Vector3& position = queryPoint.position;
  const Vector3& direction = queryPoint.direction;

  // Loop over the (currently valid) candidates and update
  for (std::size_t& index = stream.currentIndex;
       index < stream.candidates.size(); ++index) {
    // Get the candidate, and resolve the tuple
    NavigationStream::Candidate& candidate = stream.currentCandidate();
    auto& [sIntersection, portal, bTolerance, abortTarget] = candidate;
    // Get the surface from the object intersection
    const Surface* surface = sIntersection.object();
    // (re-)Intersect the surface
    auto multiIntersection = surface->intersect(
        gctx, position, direction, bTolerance, s_onSurfaceTolerance);
    // Split them into valid intersections
    for (auto& rsIntersection : multiIntersection.split()) {
      // Skip wrong index solution
      if (rsIntersection.index() != sIntersection.index()) {
        continue;
      }
      // Valid solution is either on surface or updates the distance
      if (rsIntersection.isValid()) {
        sIntersection = rsIntersection;
        return true;
      }
    }
  }
  // No candidate was reachable
  return false;
}
