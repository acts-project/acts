// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <cmath>
#include <cstddef>

std::vector<Acts::ActsScalar> Acts::detail::VerticesHelper::phiSegments(
    ActsScalar phiMin, ActsScalar phiMax,
    const std::vector<ActsScalar>& phiRefs, unsigned int minSegments) {
  // Check that the phi range is valid
  if (phiMin > phiMax) {
    throw std::invalid_argument(
        "VerticesHelper::phiSegments ... Minimum phi must be smaller than "
        "maximum phi");
  }

  // First check that no reference phi is outside the range
  for (const auto& phiRef : phiRefs) {
    if (phiRef < phiMin || phiRef > phiMax) {
      throw std::invalid_argument(
          "VerticesHelper::phiSegments ... Reference phi is outside the range "
          "of the segment");
    }
  }
  // Bail out if minSegmens is smaler 4 or not a multiple of 4
  if (minSegments < 4 || minSegments % 4 != 0) {
    throw std::invalid_argument(
        "VerticesHelper::phiSegments ...Minimum number of segments must be a "
        "multiple "
        "of 4 and at least 4");
  }
  std::vector<ActsScalar> phiSegments = {phiMin, phiMax};
  // Minimum approximation for a circle need
  for (unsigned int i = 0; i < minSegments + 1; ++i) {
    ActsScalar phiExt = -M_PI + i * 2 * M_PI / minSegments;
    if (phiExt > phiMin && phiExt < phiMax &&
        std::find_if(phiSegments.begin(), phiSegments.end(),
                     [&phiExt](ActsScalar phi) {
                       return std::abs(phi - phiExt) <
                              std::numeric_limits<ActsScalar>::epsilon();
                     }) == phiSegments.end()) {
      phiSegments.push_back(phiExt);
    }
  }
  // Add the reference phis
  for (const auto& phiRef : phiRefs) {
    if (phiRef > phiMin && phiRef < phiMax) {
      if (std::find_if(phiSegments.begin(), phiSegments.end(),
                       [&phiRef](ActsScalar phi) {
                         return std::abs(phi - phiRef) <
                                std::numeric_limits<ActsScalar>::epsilon();
                       }) == phiSegments.end()) {
        phiSegments.push_back(phiRef);
      }
    }
  }
  // Sort the phis
  std::sort(phiSegments.begin(), phiSegments.end());
  return phiSegments;
}

std::vector<Acts::Vector2> Acts::detail::VerticesHelper::ellipsoidVertices(
    ActsScalar innerRx, ActsScalar innerRy, ActsScalar outerRx,
    ActsScalar outerRy, ActsScalar avgPhi, ActsScalar halfPhi,
    unsigned int lseg) {
  // List of vertices counter-clockwise starting at smallest phi w.r.t center,
  // for both inner/outer ring/segment
  std::vector<Vector2> rvertices;  // return vertices
  std::vector<Vector2> ivertices;  // inner vertices
  std::vector<Vector2> overtices;  // outer verices

  bool innerExists = (innerRx > 0. && innerRy > 0.);
  bool closed = std::abs(halfPhi - M_PI) < s_onSurfaceTolerance;

  // The inner (if exists) and outer bow
  if (innerExists) {
    ivertices = createSegment<Vector2, Transform2>(
        {innerRx, innerRy}, avgPhi - halfPhi, avgPhi + halfPhi, {}, lseg);
  }
  overtices = createSegment<Vector2, Transform2>(
      {outerRx, outerRy}, avgPhi - halfPhi, avgPhi + halfPhi, {}, lseg);

  // We want to keep the same counter-clockwise orientation for displaying
  if (!innerExists) {
    if (!closed) {
      // Add the center case we have a sector
      rvertices.push_back(Vector2(0., 0.));
    }
    rvertices.insert(rvertices.end(), overtices.begin(), overtices.end());
  } else if (!closed) {
    rvertices.insert(rvertices.end(), overtices.begin(), overtices.end());
    rvertices.insert(rvertices.end(), ivertices.rbegin(), ivertices.rend());
  } else {
    rvertices.insert(rvertices.end(), overtices.begin(), overtices.end());
    rvertices.insert(rvertices.end(), ivertices.begin(), ivertices.end());
  }
  return rvertices;
}

std::vector<Acts::Vector2> Acts::detail::VerticesHelper::circularVertices(
    ActsScalar innerR, ActsScalar outerR, ActsScalar avgPhi, ActsScalar halfPhi,
    unsigned int lseg) {
  return ellipsoidVertices(innerR, innerR, outerR, outerR, avgPhi, halfPhi,
                           lseg);
}

bool Acts::detail::VerticesHelper::onHyperPlane(
    const std::vector<Acts::Vector3>& vertices, ActsScalar tolerance) {
  // Obvious always on one surface
  if (vertices.size() < 4) {
    return true;
  }
  // Create the hyperplane
  auto hyperPlane = Eigen::Hyperplane<ActsScalar, 3>::Through(
      vertices[0], vertices[1], vertices[2]);
  for (std::size_t ip = 3; ip < vertices.size(); ++ip) {
    if (hyperPlane.absDistance(vertices[ip]) > tolerance) {
      return false;
    }
  }
  return true;
}
