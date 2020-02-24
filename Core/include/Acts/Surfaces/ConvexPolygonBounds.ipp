// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/ThrowAssert.hpp"

std::ostream& Acts::ConvexPolygonBoundsBase::toStream(std::ostream& sl) const {
  std::vector<Vector2D> vtxs = vertices();
  sl << "Acts::ConvexPolygonBounds<" << vtxs.size() << ">: vertices: [x, y]\n";
  for (size_t i = 0; i < vtxs.size(); i++) {
    const auto& vtx = vtxs[i];
    if (i > 0) {
      sl << ",";
      sl << "\n";
    }
    sl << "[" << vtx.x() << ", " << vtx.y() << "]";
  }
  return sl;
}

template <typename coll_t>
Acts::RectangleBounds Acts::ConvexPolygonBoundsBase::makeBoundingBox(
    const coll_t& vertices) {
  Vector2D vmax, vmin;
  vmax = vertices[0];
  vmin = vertices[0];

  for (size_t i = 1; i < vertices.size(); i++) {
    vmax = vmax.cwiseMax(vertices[i]);
    vmin = vmin.cwiseMin(vertices[i]);
  }

  return {vmin, vmax};
}

std::vector<TDD_real_t> Acts::ConvexPolygonBoundsBase::valueStore() const {
  std::vector<TDD_real_t> values;
  for (const auto& vtx : vertices()) {
    values.push_back(vtx.x());
    values.push_back(vtx.y());
  }
  return values;
}

template <typename coll_t>
bool Acts::ConvexPolygonBoundsBase::convex_impl(const coll_t& vertices) {
  static_assert(std::is_same<typename coll_t::value_type, Vector2D>::value,
                "Must be collection of Vector2D");

  const size_t N = vertices.size();
  for (size_t i = 0; i < N; i++) {
    size_t j = (i + 1) % N;
    const Vector2D& a = vertices[i];
    const Vector2D& b = vertices[j];

    const Vector2D ab = b - a;
    const Vector2D normal = Vector2D(ab.y(), -ab.x()).normalized();

    bool first = true;
    bool ref;
    // loop over all other vertices
    for (size_t k = 0; k < N; k++) {
      if (k == i || k == j) {
        continue;
      }

      const Vector2D& c = vertices[k];
      double dot = normal.dot(c - a);

      if (first) {
        ref = std::signbit(dot);
        first = false;
        continue;
      }

      if (std::signbit(dot) != ref) {
        return false;
      }
    }
  }
  return true;
}

template <int N>
Acts::ConvexPolygonBounds<N>::ConvexPolygonBounds(
    const std::vector<Acts::Vector2D>& vertices)
    : m_vertices(), m_boundingBox(makeBoundingBox(vertices)) {
  throw_assert(vertices.size() == N,
               "Size and number of given vertices do not match.");
  for (size_t i = 0; i < N; i++) {
    m_vertices[i] = vertices[i];
  }
  throw_assert(convex(), "Given vertices do not form convex hull.");
}

template <int N>
Acts::ConvexPolygonBounds<N>::ConvexPolygonBounds(const vertex_array& vertices)
    : m_vertices(vertices), m_boundingBox(makeBoundingBox(vertices)) {
  throw_assert(convex(), "Given vertices do not form convex hull.");
}

template <int N>
Acts::ConvexPolygonBounds<N>* Acts::ConvexPolygonBounds<N>::clone() const {
  return new ConvexPolygonBounds<N>(*this);
}

template <int N>
Acts::SurfaceBounds::BoundsType Acts::ConvexPolygonBounds<N>::type() const {
  return SurfaceBounds::ConvexPolygon;
}

template <int N>
bool Acts::ConvexPolygonBounds<N>::inside(
    const Acts::Vector2D& lposition, const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(lposition, m_vertices);
}

template <int N>
double Acts::ConvexPolygonBounds<N>::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(lposition, m_vertices);
}

template <int N>
std::vector<Acts::Vector2D> Acts::ConvexPolygonBounds<N>::vertices(
    unsigned int /*lseg*/) const {
  return {m_vertices.begin(), m_vertices.end()};
}

template <int N>
const Acts::RectangleBounds& Acts::ConvexPolygonBounds<N>::boundingBox() const {
  return m_boundingBox;
}

template <int N>
bool Acts::ConvexPolygonBounds<N>::convex() const {
  return convex_impl(m_vertices);
}

Acts::ConvexPolygonBounds<Acts::PolygonDynamic>::ConvexPolygonBounds(
    const std::vector<Vector2D>& vertices)
    : m_vertices(vertices.begin(), vertices.end()),
      m_boundingBox(makeBoundingBox(vertices)) {}
Acts::ConvexPolygonBounds<Acts::PolygonDynamic>*
Acts::ConvexPolygonBounds<Acts::PolygonDynamic>::clone() const {
  return new Acts::ConvexPolygonBounds<Acts::PolygonDynamic>(*this);
}

Acts::SurfaceBounds::BoundsType
Acts::ConvexPolygonBounds<Acts::PolygonDynamic>::type() const {
  return SurfaceBounds::ConvexPolygon;
}

bool Acts::ConvexPolygonBounds<Acts::PolygonDynamic>::inside(
    const Acts::Vector2D& lposition, const Acts::BoundaryCheck& bcheck) const {
  return bcheck.isInside(lposition, m_vertices);
}

double Acts::ConvexPolygonBounds<Acts::PolygonDynamic>::distanceToBoundary(
    const Acts::Vector2D& lposition) const {
  return BoundaryCheck(true).distance(lposition, m_vertices);
}

std::vector<Acts::Vector2D> Acts::ConvexPolygonBounds<
    Acts::PolygonDynamic>::vertices(unsigned int /*lseg*/) const {
  return {m_vertices.begin(), m_vertices.end()};
}

const Acts::RectangleBounds&
Acts::ConvexPolygonBounds<Acts::PolygonDynamic>::boundingBox() const {
  return m_boundingBox;
}

bool Acts::ConvexPolygonBounds<Acts::PolygonDynamic>::convex() const {
  return convex_impl(m_vertices);
}
