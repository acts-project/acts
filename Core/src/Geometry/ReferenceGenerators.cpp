// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/ReferenceGenerators.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

const std::vector<Acts::Vector3> Acts::PolyhedronReferenceGenerator::references(
    const GeometryContext& gctx, const Surface& surface) const {
  // Create the return  vector
  std::vector<Vector3> rPositions;
  auto pHedron = surface.polyhedronRepresentation(gctx, nSegements);
  rPositions.insert(rPositions.end(), pHedron.vertices.begin(),
                    pHedron.vertices.end());

  // Compute the barycenter
  Vector3 bc(0., 0., 0.);
  std::ranges::for_each(rPositions, [&](const auto& p) { bc += p; });
  bc *= 1. / rPositions.size();

  // if an expansion is requested, calculate it for every vertex
  if (expansionValue != 0.0) {
    std::ranges::for_each(rPositions, [&](auto& p) {
      p += expansionValue * Vector3(p - bc).normalized();
    });
  }

  // Add the barycenter if configured
  if (addBarycenter) {
    rPositions.push_back(bc);
  }
  return rPositions;
}

const std::vector<Acts::Vector3> Acts::ProjectedReferenceGenerator::references(
    const GeometryContext& gctx, const Surface& surface) const {
  if (referenceSurface == nullptr) {
    throw std::invalid_argument(
        "ProjectedReferenceGenerator: reference surface is nullptr.");
  }

  // Create the test polyhedron
  std::vector<Vector3> tPositions;
  auto pHedron = surface.polyhedronRepresentation(gctx, nSegements);
  tPositions.insert(tPositions.end(), pHedron.vertices.begin(),
                    pHedron.vertices.end());

  // Create intersected position
  std::vector<Vector3> rPositions;
  Vector3 rCog;  // reference center of gravity
  // Loop over luminous regious points and project
  for (const auto& lp : luminousRegion) {
    for (const auto& tp : tPositions) {
      // Create the ray from luminous point to test point
      Vector3 rayDirection = (tp - lp).normalized();
      auto refMultiIntersections =
          referenceSurface->intersect(gctx, lp, rayDirection);
      // Take the closest intersection point in forward direction
      rPositions.push_back(refMultiIntersections.closestForward().position());
      rCog += rPositions.back();
    }
  }

  // Normalize center of gravity
  rCog /= rPositions.size();

  // if an expansion is requested, calculate it for every vertex
  if (expansionValue != 0.0) {
    std::ranges::for_each(rPositions, [&](auto& p) {
      p += expansionValue * Vector3(p - rCog).normalized();
    });
  }

  return rPositions;
}
