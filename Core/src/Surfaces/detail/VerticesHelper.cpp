// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

namespace Acts {

std::vector<double> detail::VerticesHelper::phiSegments(
    double phiMin, double phiMax, const std::vector<double>& phiRefs,
    unsigned int quarterSegments) {
  // Check that the phi range is valid
  if (phiMin > phiMax) {
    throw std::invalid_argument(
        "VerticesHelper::phiSegments ... Minimum phi must be smaller than "
        "maximum phi");
  }

  // First check that no reference phi is outside the range
  for (double phiRef : phiRefs) {
    if (phiRef < phiMin || phiRef > phiMax) {
      throw std::invalid_argument(
          "VerticesHelper::phiSegments ... Reference phi is outside the range "
          "of the segment");
    }
  }
  if (quarterSegments == 0u) {
    throw std::invalid_argument(
        "VerticesHelper::phiSegments ... Number of segments must be larger "
        "than 0.");
  }
  std::vector<double> phiSegments = {phiMin, phiMax};
  // Minimum approximation for a circle need
  // - if the circle is closed the last point is given twice
  for (unsigned int i = 0; i < 4 * quarterSegments + 1; ++i) {
    double phiExt =
        -std::numbers::pi + i * 2 * std::numbers::pi / (4 * quarterSegments);
    if (phiExt > phiMin && phiExt < phiMax &&
        std::ranges::none_of(phiSegments, [&phiExt](double phi) {
          return std::abs(phi - phiExt) <
                 std::numeric_limits<double>::epsilon();
        })) {
      phiSegments.push_back(phiExt);
    }
  }
  // Add the reference phis
  for (const auto& phiRef : phiRefs) {
    if (phiRef > phiMin && phiRef < phiMax) {
      if (std::ranges::none_of(phiSegments, [&phiRef](double phi) {
            return std::abs(phi - phiRef) <
                   std::numeric_limits<double>::epsilon();
          })) {
        phiSegments.push_back(phiRef);
      }
    }
  }

  // Sort the phis
  std::ranges::sort(phiSegments);
  return phiSegments;
}

std::vector<Vector2> detail::VerticesHelper::ellipsoidVertices(
    double innerRx, double innerRy, double outerRx, double outerRy,
    double avgPhi, double halfPhi, unsigned int quarterSegments) {
  // List of vertices counter-clockwise starting at smallest phi w.r.t center,
  // for both inner/outer ring/segment
  std::vector<Vector2> rvertices;  // return vertices
  std::vector<Vector2> ivertices;  // inner vertices
  std::vector<Vector2> overtices;  // outer verices

  bool innerExists = (innerRx > 0. && innerRy > 0.);
  bool closed = std::abs(halfPhi - std::numbers::pi) < s_onSurfaceTolerance;

  std::vector<double> refPhi = {};
  if (avgPhi != 0.) {
    refPhi.push_back(avgPhi);
  }

  // The inner (if exists) and outer bow
  if (innerExists) {
    ivertices = segmentVertices<Vector2, Transform2>(
        {innerRx, innerRy}, avgPhi - halfPhi, avgPhi + halfPhi, refPhi,
        quarterSegments);
  }
  overtices = segmentVertices<Vector2, Transform2>(
      {outerRx, outerRy}, avgPhi - halfPhi, avgPhi + halfPhi, refPhi,
      quarterSegments);

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

std::vector<Vector2> detail::VerticesHelper::circularVertices(
    double innerR, double outerR, double avgPhi, double halfPhi,
    unsigned int quarterSegments) {
  return ellipsoidVertices(innerR, innerR, outerR, outerR, avgPhi, halfPhi,
                           quarterSegments);
}

bool detail::VerticesHelper::onHyperPlane(const std::vector<Vector3>& vertices,
                                          double tolerance) {
  // Obvious always on one surface
  if (vertices.size() < 4) {
    return true;
  }
  // Create the hyperplane
  auto hyperPlane = Eigen::Hyperplane<double, 3>::Through(
      vertices[0], vertices[1], vertices[2]);
  for (std::size_t ip = 3; ip < vertices.size(); ++ip) {
    if (hyperPlane.absDistance(vertices[ip]) > tolerance) {
      return false;
    }
  }
  return true;
}

Vector2 detail::VerticesHelper::computeClosestPointOnPolygon(
    const Vector2& point, std::span<const Vector2> vertices,
    const SquareMatrix2& metric) {
  auto squaredNorm = [&](const Vector2& x) {
    return (x.transpose() * metric * x).value();
  };

  // calculate the closest position on the segment between `ll0` and `ll1` to
  // the point as measured by the metric induced by the metric matrix
  auto closestOnSegment = [&](auto&& ll0, auto&& ll1) {
    // normal vector and position of the closest point along the normal
    auto n = ll1 - ll0;
    auto n_transformed = metric * n;
    auto f = n.dot(n_transformed);
    auto u = std::isnormal(f)
                 ? (point - ll0).dot(n_transformed) / f
                 : 0.5;  // ll0 and ll1 are so close it doesn't matter
    // u must be in [0, 1] to still be on the polygon segment
    return ll0 + std::clamp(u, 0.0, 1.0) * n;
  };

  auto iv = std::begin(vertices);
  Vector2 l0 = *iv;
  Vector2 l1 = *(++iv);
  Vector2 closest = closestOnSegment(l0, l1);
  auto closestDist = squaredNorm(closest - point);
  // Calculate the closest point on other connecting lines and compare distances
  for (++iv; iv != std::end(vertices); ++iv) {
    l0 = l1;
    l1 = *iv;
    Vector2 current = closestOnSegment(l0, l1);
    auto currentDist = squaredNorm(current - point);
    if (currentDist < closestDist) {
      closest = current;
      closestDist = currentDist;
    }
  }
  // final edge from last vertex back to the first vertex
  Vector2 last = closestOnSegment(l1, *std::begin(vertices));
  if (squaredNorm(last - point) < closestDist) {
    closest = last;
  }
  return closest;
}

Vector2 detail::VerticesHelper::computeEuclideanClosestPointOnRectangle(
    const Vector2& point, const Vector2& lowerLeft, const Vector2& upperRight) {
  /*
   *
   *        |                 |
   *   IV   |       V         | I
   *        |                 |
   *  ------------------------------
   *        |                 |
   *        |                 |
   *   VIII |     INSIDE      | VI
   *        |                 |
   *        |                 |
   *  ------------------------------
   *        |                 |
   *   III  |      VII        | II
   *        |                 |
   *
   */

  double l0 = point[0];
  double l1 = point[1];
  double loc0Min = lowerLeft[0];
  double loc0Max = upperRight[0];
  double loc1Min = lowerLeft[1];
  double loc1Max = upperRight[1];

  // check if inside
  if (loc0Min <= l0 && l0 < loc0Max && loc1Min <= l1 && l1 < loc1Max) {
    // INSIDE
    double dist = std::abs(loc0Max - l0);
    Vector2 cls(loc0Max, l1);

    double test = std::abs(loc0Min - l0);
    if (test <= dist) {
      dist = test;
      cls = {loc0Min, l1};
    }

    test = std::abs(loc1Max - l1);
    if (test <= dist) {
      dist = test;
      cls = {l0, loc1Max};
    }

    test = std::abs(loc1Min - l1);
    if (test <= dist) {
      return {l0, loc1Min};
    }
    return cls;
  } else {
    // OUTSIDE, check sectors
    if (l0 > loc0Max) {
      if (l1 > loc1Max) {  // I
        return {loc0Max, loc1Max};
      } else if (l1 <= loc1Min) {  // II
        return {loc0Max, loc1Min};
      } else {  // VI
        return {loc0Max, l1};
      }
    } else if (l0 < loc0Min) {
      if (l1 > loc1Max) {  // IV
        return {loc0Min, loc1Max};
      } else if (l1 <= loc1Min) {  // III
        return {loc0Min, loc1Min};
      } else {  // VIII
        return {loc0Min, l1};
      }
    } else {
      if (l1 > loc1Max) {  // V
        return {l0, loc1Max};
      } else {  // l1 <= loc1Min # VII
        return {l0, loc1Min};
      }
      // third case not necessary, see INSIDE above
    }
  }
}

Vector2 detail::VerticesHelper::computeClosestPointOnAlignedBox(
    const Vector2& lowerLeft, const Vector2& upperRight, const Vector2& point,
    const SquareMatrix2& metric) {
  Vector2 closestPoint;

  if (metric.isIdentity()) {
    closestPoint =
        detail::VerticesHelper::computeEuclideanClosestPointOnRectangle(
            point, lowerLeft, upperRight);
  } else {
    // TODO there might be a more optimal way to compute the closest point to a
    // box with metric

    std::array<Vector2, 4> vertices = {{lowerLeft,
                                        {upperRight[0], lowerLeft[1]},
                                        upperRight,
                                        {lowerLeft[0], upperRight[1]}}};

    closestPoint = detail::VerticesHelper::computeClosestPointOnPolygon(
        point, vertices, metric);
  }

  return closestPoint;
}

}  // namespace Acts
