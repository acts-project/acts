// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <utility>
#include <vector>
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

namespace detail {

/// @brief Helper method for set of vertices for polyhedrons
/// drawing and inside/outside checks
namespace VerticesHelper {

/// A method that inserts the cartesian extrema points and segments
/// a curved segment into sub segments
///
/// @param phiMin the minimum Phi of the bounds object
/// @param phiMax the maximum Phi of the bounds object
/// @param phiRef is a vector of reference phi values to be included as well
/// @param phiTolerance is the tolerance for reference phi insertion
///
/// @return a vector
static std::vector<double> phiSegments(double phiMin = -M_PI,
                                       double phiMax = M_PI,
                                       std::vector<double> phiRefs = {},
                                       double phiTolerance = 1e-6) {
  // This is to ensure that the extrema are built regardless of number
  // of segments
  std::vector<double> phiSegments;
  std::vector<double> quarters = {-M_PI, -0.5 * M_PI, 0., 0.5 * M_PI, M_PI};
  // It does not cover the full azimuth
  if (phiMin != -M_PI or phiMax != M_PI) {
    phiSegments.push_back(phiMin);
    for (unsigned int iq = 1; iq < 4; ++iq) {
      if (phiMin < quarters[iq] and phiMax > quarters[iq]) {
        phiSegments.push_back(quarters[iq]);
      }
    }
    phiSegments.push_back(phiMax);
  } else {
    phiSegments = quarters;
  }
  // Insert the reference phis if
  if (not phiRefs.empty()) {
    for (const auto& phiRef : phiRefs) {
      // Trying to find the right patch
      auto match = std::find_if(
          phiSegments.begin(), phiSegments.end(), [&](double phiSeg) {
            return std::abs(phiSeg - phiRef) < phiTolerance;
          });
      if (match == phiSegments.end()) {
        phiSegments.push_back(phiRef);
      }
    }
    std::sort(phiSegments.begin(), phiSegments.end());
  }
  return phiSegments;
}

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
                   std::pair<double, double> rxy, double phi1, double phi2,
                   unsigned int lseg, int addon = 0,
                   const vertex_t& offset = vertex_t::Zero(),
                   const transform_t& transform = transform_t::Identity()) {
  // Calculate the number of segments - 1 is the minimum
  unsigned int segs = std::abs(phi2 - phi1) / (2 * M_PI) * lseg;
  segs = segs > 0 ? segs : 1;
  double phistep = (phi2 - phi1) / segs;
  // Create the segments
  for (unsigned int iphi = 0; iphi < segs + addon; ++iphi) {
    double phi = phi1 + iphi * phistep;
    vertex_t vertex = vertex_t::Zero();
    vertex(0) = rxy.first * std::cos(phi);
    vertex(1) = rxy.second * std::sin(phi);
    vertex = vertex + offset;
    vertices.push_back(transform * vertex);
  }
}

/// Vertices on an ellipse-like bound object
///
/// @param innerRx The radius of the inner ellipse (in x), 0 if sector
/// @param innerRy The radius of the inner ellipse (in y), 0 if sector
/// @param outerRx The radius of the outer ellipse (in x)
/// @param outerRy The radius of the outer ellipse (in y)
/// @param avgPhi The phi direction of the center if sector
/// @param halfPhi The half phi sector if sector
/// @param lseg The number of segments for for a full 2*pi segment
///
/// @return a vector of 2d-vectors
static std::vector<Vector2D> ellispoidVertices(double innerRx, double innerRy,
                                               double outerRx, double outerRy,
                                               double avgPhi = 0.,
                                               double halfPhi = M_PI,
                                               unsigned int lseg = 1) {
  // List of vertices counter-clockwise starting at smallest phi w.r.t center,
  // for both inner/outer ring/segment
  std::vector<Acts::Vector2D> rvertices;  // return vertices
  std::vector<Acts::Vector2D> ivertices;  // inner vertices
  std::vector<Acts::Vector2D> overtices;  // outer verices

  bool innerExists = (innerRx > 0. and innerRy > 0.);
  bool closed = std::abs(halfPhi - M_PI) < s_onSurfaceTolerance;

  // Get the phi segments from the helper method
  auto phiSegs = detail::VerticesHelper::phiSegments(
      avgPhi - halfPhi, avgPhi + halfPhi, {avgPhi});

  // The inner (if exists) and outer bow
  for (unsigned int iseg = 0; iseg < phiSegs.size() - 1; ++iseg) {
    int addon = (iseg == phiSegs.size() - 2 and not closed) ? 1 : 0;
    if (innerExists) {
      detail::VerticesHelper::createSegment<Vector2D, Eigen::Affine2d>(
          ivertices, {innerRx, innerRy}, phiSegs[iseg], phiSegs[iseg + 1], lseg,
          addon);
    }
    detail::VerticesHelper::createSegment<Vector2D, Eigen::Affine2d>(
        overtices, {outerRx, outerRy}, phiSegs[iseg], phiSegs[iseg + 1], lseg,
        addon);
  }

  // We want to keep the same counter-clockwise orientation for displaying
  if (not innerExists) {
    if (not closed) {
      // Add the center case we have a sector
      rvertices.push_back(Vector2D(0., 0.));
    }
    rvertices.insert(rvertices.end(), overtices.begin(), overtices.end());
  } else if (not closed) {
    rvertices.insert(rvertices.end(), overtices.begin(), overtices.end());
    rvertices.insert(rvertices.end(), ivertices.rbegin(), ivertices.rend());
  } else {
    rvertices.insert(rvertices.end(), overtices.begin(), overtices.end());
    rvertices.insert(rvertices.end(), ivertices.begin(), ivertices.end());
  }
  return rvertices;
}

/// Vertices on an disc/wheel-like bound object
///
/// @param innerR The radius of the inner circle (sector)
/// @param outerR The radius of the outer circle (sector)
/// @param avgPhi The phi direction of the center if sector
/// @param halfPhi The half phi sector if sector
/// @param lseg The number of segments for for a full 2*pi segment
///
/// @return a vector of 2d-vectors
static std::vector<Vector2D> circularVertices(double innerR, double outerR,
                                              double avgPhi = 0.,
                                              double halfPhi = M_PI,
                                              unsigned int lseg = 1) {
  return ellispoidVertices(innerR, innerR, outerR, outerR, avgPhi, halfPhi,
                           lseg);
}

/// Check if the point is inside the polygon w/o any tolerances
///
/// @tparam vertex_container_t is an iterable container
///
/// @param point is the Vector2DType to check
/// @param vertices Forward iterable container of convex polygon vertices.
///                 Calling `std::begin`/ `std::end` on the container must
///                 return an iterator where `*it` must be convertible to
///                 an `Vector2DType`.
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

/// Check if the point is inside the rectangle
///
/// @tparam vertex_t is vector with [0],[1] access
///
/// @param point is the Vector2DType to check
/// @param vertices Forward iterable container of convex polygon vertices.
///                 Calling `std::begin`/ `std::end` on the container must
///                 return an iterator where `*it` must be convertible to
///                 an `Vector2DType`.
/// @return bool for inside/outside
template <typename vertex_t>
bool isInsideRectangle(const vertex_t& point, const vertex_t& lowerLeft,
                       const vertex_t& upperRight) {
  return (lowerLeft[0] <= point[0]) && (point[0] < upperRight[0]) &&
         (lowerLeft[1] <= point[1]) && (point[1] < upperRight[1]);
}

};  // namespace VerticesHelper

}  // namespace detail

}  // namespace Acts