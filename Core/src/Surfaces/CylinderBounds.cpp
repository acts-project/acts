// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/CylinderBounds.hpp"

#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Utilities/detail/OstreamStateGuard.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <utility>

namespace Acts {

using VectorHelpers::perp;
using VectorHelpers::phi;

std::vector<double> CylinderBounds::values() const {
  return {m_values.begin(), m_values.end()};
}

Vector2 CylinderBounds::shifted(const Vector2& lposition) const {
  return {detail::radian_sym((lposition[0] / get(eR)) - get(eAveragePhi)),
          lposition[1]};
}

bool CylinderBounds::inside(const Vector2& lposition) const {
  double halfLengthZ = get(eHalfLengthZ);
  double halfPhi = get(eHalfPhiSector);
  return detail::VerticesHelper::isInsideRectangle(
      shifted(lposition), Vector2(-halfPhi, -halfLengthZ),
      Vector2(halfPhi, halfLengthZ));
}

Vector2 CylinderBounds::closestPoint(const Vector2& lposition,
                                     const SquareMatrix2& metric) const {
  double halfLengthZ = get(eHalfLengthZ);
  double radius = get(eR);

  Vector2 vertices[] = {{-radius, -halfLengthZ},
                        {radius, -halfLengthZ},
                        {radius, halfLengthZ},
                        {-radius, halfLengthZ}};

  return detail::VerticesHelper::computeClosestPointOnPolygon(lposition,
                                                              vertices, metric);
}

Vector2 CylinderBounds::center() const {
  // For cylinder bounds in local coordinates (rphi, z),
  // centroid is at (averagePhi, 0) since z extends symmetrically
  return Vector2(get(eAveragePhi), 0.0);
}

std::ostream& CylinderBounds::toStream(std::ostream& sl) const {
  detail::OstreamStateGuard guard{sl};
  sl << std::fixed << std::setprecision(7);
  sl << "Acts::CylinderBounds: (radius, halfLengthZ, halfPhiSector, "
        "averagePhi) = ";
  sl << "(" << get(eR) << ", " << get(eHalfLengthZ) << ", ";
  sl << get(eHalfPhiSector) << ", " << get(eAveragePhi) << ")";
  return sl;
}

std::vector<Vector3> CylinderBounds::circleVertices(
    const Transform3 transform, unsigned int quarterSegments) const {
  std::vector<Vector3> vertices;

  double avgPhi = get(eAveragePhi);
  double halfPhi = get(eHalfPhiSector);

  std::vector<double> phiRef = {};
  if (bool fullCylinder = coversFullAzimuth(); fullCylinder) {
    phiRef = {avgPhi};
  }

  // Write the two bows/circles on either side
  std::vector<int> sides = {-1, 1};
  for (auto& side : sides) {
    // Helper method to create the segment
    auto svertices = detail::VerticesHelper::segmentVertices(
        {get(eR), get(eR)}, avgPhi - halfPhi, avgPhi + halfPhi, phiRef,
        quarterSegments, Vector3(0., 0., side * get(eHalfLengthZ)), transform);
    vertices.insert(vertices.end(), svertices.begin(), svertices.end());
  }

  return vertices;
}

void CylinderBounds::checkConsistency() noexcept(false) {
  if (get(eR) <= 0.) {
    throw std::invalid_argument(
        "CylinderBounds: invalid radial setup: radius is negative");
  }
  if (get(eHalfLengthZ) <= 0.) {
    throw std::invalid_argument(
        "CylinderBounds: invalid length setup: half length is negative");
  }
  if (get(eHalfPhiSector) <= 0. || get(eHalfPhiSector) > std::numbers::pi) {
    throw std::invalid_argument("CylinderBounds: invalid phi sector setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi)) &&
      std::abs(std::abs(get(eAveragePhi)) - std::numbers::pi) > s_epsilon) {
    throw std::invalid_argument("CylinderBounds: invalid phi positioning.");
  }
}

}  // namespace Acts
