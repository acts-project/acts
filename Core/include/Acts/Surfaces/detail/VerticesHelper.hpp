// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

/// Helper methods for polyhedron vertices drawing and inside/outside checks.
namespace Acts::detail::VerticesHelper {

/// A method that inserts the cartesian extrema points and segments
/// a curved segment into sub segments
///
/// @param phiMin the minimum Phi of the bounds object
/// @param phiMax the maximum Phi of the bounds object
/// @param phiRef is a vector of reference phi values to be included as well
/// @param phiTolerance is the tolerance for reference phi insertion
/// @return a vector
std::vector<ActsScalar> phiSegments(ActsScalar phiMin = -M_PI,
                                    ActsScalar phiMax = M_PI,
                                    const std::vector<ActsScalar>& phiRefs = {},
                                    ActsScalar phiTolerance = 1e-6);

/// Helper method to create a regular 2 or 3 D segment
///  between two phi values
///
/// @tparam vertex_t Type of vertex to be applied
/// @tparam transform_t Optional transform
///
/// @param vertices [in,out] The 3D vertices to be filled
/// @param rxy The radius description if first +/= second: ellipse
/// @param phi1 The first phi value
/// @param phi2 The second phi value
/// @param lseg The number of segments for full 2*PI
/// @param addon The additional segments to be built
/// @param offset The out of plane offset position of the bow
/// @param transform The transform applied (optional)
template <typename vertex_t, typename transform_t>
void createSegment(std::vector<vertex_t>& vertices,
                   std::pair<ActsScalar, ActsScalar> rxy, ActsScalar phi1,
                   ActsScalar phi2, unsigned int lseg, int addon = 0,
                   const vertex_t& offset = vertex_t::Zero(),
                   const transform_t& transform = transform_t::Identity()) {
  // Calculate the number of segments - 1 is the minimum
  unsigned int segs =
      static_cast<unsigned int>(std::abs(phi2 - phi1) / (2 * M_PI) * lseg);
  segs = segs > 0 ? segs : 1;
  ActsScalar phistep = (phi2 - phi1) / segs;
  // Create the segments
  for (unsigned int iphi = 0; iphi < segs + addon; ++iphi) {
    ActsScalar phi = phi1 + iphi * phistep;
    vertex_t vertex = vertex_t::Zero();
    vertex(0) = rxy.first * std::cos(phi);
    vertex(1) = rxy.second * std::sin(phi);

    vertex = vertex + offset;
    vertices.push_back(transform * vertex);
  }
}

/// Construct vertices on an ellipse-like bound object.
///
/// @param innerRx The radius of the inner ellipse (in x), 0 if sector
/// @param innerRy The radius of the inner ellipse (in y), 0 if sector
/// @param outerRx The radius of the outer ellipse (in x)
/// @param outerRy The radius of the outer ellipse (in y)
/// @param avgPhi The phi direction of the center if sector
/// @param halfPhi The half phi sector if sector
/// @param lseg The number of segments for for a full 2*pi segment
/// @return a vector of 2d-vectors
std::vector<Vector2> ellipsoidVertices(ActsScalar innerRx, ActsScalar innerRy,
                                       ActsScalar outerRx, ActsScalar outerRy,
                                       ActsScalar avgPhi = 0.,
                                       ActsScalar halfPhi = M_PI,
                                       unsigned int lseg = 1);

/// Construct vertices on an disc/wheel-like bound object.
///
/// @param innerR The radius of the inner circle (sector)
/// @param outerR The radius of the outer circle (sector)
/// @param avgPhi The phi direction of the center if sector
/// @param halfPhi The half phi sector if sector
/// @param lseg The number of segments for for a full 2*pi segment
/// @return a vector of 2d-vectors
std::vector<Vector2> circularVertices(ActsScalar innerR, ActsScalar outerR,
                                      ActsScalar avgPhi = 0.,
                                      ActsScalar halfPhi = M_PI,
                                      unsigned int lseg = 1);
/// Check if the point is inside the polygon w/o any tolerances.
///
/// @tparam vertex_container_t is an iterable container
///
/// @param point is the Vector2Type to check
/// @param vertices Forward iterable container of convex polygon vertices.
///                 Calling `std::begin`/ `std::end` on the container must
///                 return an iterator where `*it` must be convertible to
///                 an `Vector2Type`.
/// @return bool for inside/outside
template <typename vertex_t, typename vertex_container_t>
bool isInsidePolygon(const vertex_t& point,
                     const vertex_container_t& vertices) {
  // when we move along the edges of a convex polygon, a point on the inside of
  // the polygon will always appear on the same side of each edge.
  // a point on the outside will switch sides at least once.

  // returns which side of the connecting line between `ll0` and `ll1` the point
  // `p` is on. computes the sign of the z-component of the cross-product
  // between the line normal vector and the vector from `ll0` to `p`.
  auto lineSide = [&](auto&& ll0, auto&& ll1) {
    auto normal = ll1 - ll0;
    auto delta = point - ll0;
    return std::signbit((normal[0] * delta[1]) - (normal[1] * delta[0]));
  };

  auto iv = std::begin(vertices);
  auto l0 = *iv;
  auto l1 = *(++iv);
  // use vertex0 to vertex1 to define reference sign and compare w/ all edges
  auto reference = lineSide(l0, l1);
  for (++iv; iv != std::end(vertices); ++iv) {
    l0 = l1;
    l1 = *iv;
    if (lineSide(l0, l1) != reference) {
      return false;
    }
  }
  // manual check for last edge from last vertex back to the first vertex
  if (lineSide(l1, *std::begin(vertices)) != reference) {
    return false;
  }
  // point was always on the same side. point must be inside.
  return true;
}

/// Check if the point is inside the rectangle.
///
/// @tparam vertex_t is vector with [0],[1] access
///
/// @param point is the Vector2Type to check
/// @param vertices Forward iterable container of convex polygon vertices.
///                 Calling `std::begin`/ `std::end` on the container must
///                 return an iterator where `*it` must be convertible to
///                 an `Vector2Type`.
/// @return bool for inside/outside
template <typename vertex_t>
bool isInsideRectangle(const vertex_t& point, const vertex_t& lowerLeft,
                       const vertex_t& upperRight) {
  return (lowerLeft[0] <= point[0]) && (point[0] < upperRight[0]) &&
         (lowerLeft[1] <= point[1]) && (point[1] < upperRight[1]);
}

/// This method checks if a cloud of points are on 2D hyper-plane in 3D space.
///
/// @param vertices The list of vertices to test
/// @param tolerance The allowed out of plane tolerance
/// @return boolean to indicate if all points are inside/outside
bool onHyperPlane(const std::vector<Vector3>& vertices,
                  ActsScalar tolerance = s_onSurfaceTolerance);

}  // namespace Acts::detail::VerticesHelper
