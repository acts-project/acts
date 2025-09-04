// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/ConvexPolygonBounds.hpp"

#include "Acts/Surfaces/detail/VerticesHelper.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

namespace Acts {

template <int N>
  requires isValidConvexPolygonSize<N>
ConvexPolygonBounds<N>::ConvexPolygonBounds(
    std::span<const Vector2> vertices) noexcept(false) {
  throw_assert(vertices.size() == N,
               "Size and number of given vertices do not match.");

  std::ranges::copy(vertices, m_vertices.begin());

  checkConsistency(m_vertices);
  makeBoundingBox(m_vertices);
  calculateCenter(m_vertices);
}

template <int N>
  requires isValidConvexPolygonSize<N>
ConvexPolygonBounds<N>::ConvexPolygonBounds(
    std::span<const double> values) noexcept(false) {
  throw_assert(values.size() == 2 * N,
               "Size and number of given values do not match.");
  for (std::size_t i = 0; i < N; i++) {
    m_vertices[i] = Vector2(values[2 * i], values[2 * i + 1]);
  }
  checkConsistency(m_vertices);
  makeBoundingBox(m_vertices);
  calculateCenter(m_vertices);
}

template <int N>
  requires isValidConvexPolygonSize<N>
bool ConvexPolygonBounds<N>::inside(const Vector2& lposition) const {
  return detail::VerticesHelper::isInsidePolygon(lposition, m_vertices);
}

template <int N>
  requires isValidConvexPolygonSize<N>
Vector2 ConvexPolygonBounds<N>::closestPoint(
    const Vector2& lposition, const SquareMatrix2& metric) const {
  return detail::VerticesHelper::computeClosestPointOnPolygon(
      lposition, m_vertices, metric);
}

template <int N>
  requires isValidConvexPolygonSize<N>
std::vector<Vector2> ConvexPolygonBounds<N>::vertices(
    unsigned int /*ignoredSegments*/) const {
  return {m_vertices.begin(), m_vertices.end()};
}

}  // namespace Acts
