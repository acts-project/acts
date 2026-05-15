// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/ConvexPolygonBounds.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <ostream>
#include <stdexcept>

namespace Acts {

std::ostream& ConvexPolygonBoundsBase::toStream(std::ostream& sl) const {
  std::vector<Vector2> vtxs = vertices();
  sl << "Acts::ConvexPolygonBounds<" << vtxs.size() << ">: vertices: [x, y]\n";
  for (std::size_t i = 0; i < vtxs.size(); i++) {
    const auto& vtx = vtxs[i];
    if (i > 0) {
      sl << ",";
      sl << "\n";
    }
    sl << "[" << vtx.x() << ", " << vtx.y() << "]";
  }
  return sl;
}

std::vector<double> ConvexPolygonBoundsBase::values() const {
  std::vector<double> values;
  for (const auto& vtx : vertices()) {
    values.push_back(vtx.x());
    values.push_back(vtx.y());
  }
  return values;
}

void ConvexPolygonBoundsBase::calculateCenter(
    std::span<const Vector2> vertices) {
  Vector2 sum = Vector2::Zero();
  for (const auto& vertex : vertices) {
    sum += vertex;
  }
  m_center = sum / static_cast<double>(vertices.size());
}

Vector2 ConvexPolygonBoundsBase::center() const {
  return m_center;
}

const RectangleBounds& ConvexPolygonBoundsBase::boundingBox() const {
  return m_boundingBox;
}

void ConvexPolygonBoundsBase::makeBoundingBox(
    std::span<const Vector2> vertices) {
  Vector2 vmax, vmin;
  vmax = vertices[0];
  vmin = vertices[0];

  for (std::size_t i = 1; i < vertices.size(); i++) {
    vmax = vmax.cwiseMax(vertices[i]);
    vmin = vmin.cwiseMin(vertices[i]);
  }

  m_boundingBox = {vmin, vmax};
}

void ConvexPolygonBoundsBase::checkConsistency(
    std::span<const Vector2> vertices) noexcept(false) {
  const std::size_t N = vertices.size();
  for (std::size_t i = 0; i < N; i++) {
    std::size_t j = (i + 1) % N;
    const Vector2& a = vertices[i];
    const Vector2& b = vertices[j];

    const Vector2 ab = b - a;
    const Vector2 normal = Vector2(ab.y(), -ab.x()).normalized();

    bool first = true;
    bool ref = false;
    // loop over all other vertices
    for (std::size_t k = 0; k < N; k++) {
      if (k == i || k == j) {
        continue;
      }

      const Vector2& c = vertices[k];
      double dot = normal.dot(c - a);

      if (first) {
        ref = std::signbit(dot);
        first = false;
        continue;
      }

      if (std::signbit(dot) != ref) {
        throw std::logic_error(
            "ConvexPolygon: Given vertices do not form convex hull");
      }
    }
  }
}

ConvexPolygonBounds<PolygonDynamic>::ConvexPolygonBounds(
    std::span<const Vector2> vertices)
    : m_vertices(vertices.begin(), vertices.end()) {
  if (vertices.size() < 3) {
    throw std::invalid_argument(
        "ConvexPolygonBounds: At least 3 vertices are required.");
  }
  checkConsistency(vertices);
  calculateCenter(vertices);
  makeBoundingBox(vertices);
}

ConvexPolygonBounds<PolygonDynamic>::ConvexPolygonBounds(
    const std::vector<Vector2>& vertices)
    : ConvexPolygonBounds{std::span<const Vector2>{vertices}} {}

bool ConvexPolygonBounds<PolygonDynamic>::inside(
    const Vector2& lposition) const {
  return detail::VerticesHelper::isInsidePolygon(lposition, m_vertices);
}

Vector2 ConvexPolygonBounds<PolygonDynamic>::closestPoint(
    const Vector2& lposition, const SquareMatrix2& metric) const {
  return detail::VerticesHelper::computeClosestPointOnPolygon(
      lposition, std::span<const Vector2>(m_vertices.data(), m_vertices.size()),
      metric);
}

std::vector<Vector2> ConvexPolygonBounds<PolygonDynamic>::vertices(
    unsigned int /*lseg*/) const {
  return {m_vertices.begin(), m_vertices.end()};
}

}  // namespace Acts
