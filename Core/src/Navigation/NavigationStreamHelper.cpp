// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/NavigationStreamHelper.hpp"

#include "Acts/Surfaces/Surface.hpp"

bool Acts::NavigationStreamHelper::processStream(NavigationStream& stream,
                                                const GeometryContext& gctx,
                                                const Vector3& position,
                                                const Vector3& direction,
                                                bool initial,
                                                BoundaryTolerance cTolerance) {
  // Loop over the (currently valid) candidates
  //
  // Except for the initial update, the loop is stopped at the first reachable
  // surface.
  //
  for (size_t& index = stream.currentIndex;  index < stream.candidates.size(); ++index) {

    std::cout << "Processing candidate " << index << std::endl;
    // Get the candidate, and resolve the tuple
    NavigationStream::Candidate& candidate = stream.currentCandidate();
    auto& [sIntersection, portal, bTolerance] = candidate;
    // Take the candidate Tolerance in case it is an intial update and not a portal
    const BoundaryTolerance& tolerance = (portal || !initial) ? bTolerance : cTolerance;

    // Get the surface from the object intersection
    const Surface* surface = sIntersection.object();
    // (re-)Intersect the surface
    auto multiIntersection = surface->intersect(
        gctx, position, direction, tolerance, s_onSurfaceTolerance);

    // Split them into valid intersections
    for (auto& rsIntersection : multiIntersection.split()) {
      // Skip negative solutions for inital update, overstepping can not happen
      if (initial && rsIntersection.pathLength() < 0.) {
        continue;
      }

      // Skip wrong index solution for non-inital updates
      if (!initial && rsIntersection.index() != sIntersection.index()) {
        continue;
      }

      // Valid solution is either on surface or updates the distance
      if (rsIntersection.isValid()) {
        std::cout << " is valid " << std::endl;
        // Valid intersection, we assume ordering, update
        sIntersection = rsIntersection;
        if (!initial) {
          // We are done with the current candidate
          return true;
        }
      }
    }
  }

  // In the initial update case, we need to sort and estimate the range
  std::sort(
      stream.candidates.begin(), stream.candidates.end(),
      [](const NavigationStream::Candidate& a,
         const NavigationStream::Candidate& b) {
        const auto& [aIntersection, aPortal, aTolerance] = a;
        const auto& [bIntersection, bPortal, bTolerance] = b;
        return (aIntersection.pathLength() < bIntersection.pathLength()) &&
               aIntersection.isValid();
      });

  // The we find the first invalid candidate
  auto firstInvalid =
      std::find_if(stream.candidates.begin(), stream.candidates.end(),
                   [](const NavigationStream::Candidate& a) {
                     const auto& [aIntersection, aPortal, aTolerance] = a;
                     return !aIntersection.isValid();
                   });

  // Set the range and initialize
  stream.candidates.resize(std::distance(stream.candidates.begin(), firstInvalid));
  if (stream.candidates.empty()) {
    return false;
  }
  stream.currentIndex = 0;
  return true;
}
