// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/detail/VerticesHelper.hpp"

#include <limits>
#include <utility>
#include <variant>

namespace Acts {

template <typename T>
bool BoundaryTolerance::holdsVariant() const {
  return std::holds_alternative<T>(m_variant);
}

template <typename T>
const T& BoundaryTolerance::getVariant() const {
  return std::get<T>(m_variant);
}

template <typename T>
const T* BoundaryTolerance::getVariantPtr() const {
  return holdsVariant<T>() ? &getVariant<T>() : nullptr;
}

template <typename Vector2Container>
PolygonBoundaryCheck<Vector2Container>::PolygonBoundaryCheck(
    const Vector2Container& vertices, BoundaryTolerance tolerance)
    : m_vertices{vertices}, m_tolerance{std::move(tolerance)} {}

template <typename Vector2Container>
const Vector2Container& PolygonBoundaryCheck<Vector2Container>::vertices()
    const {
  return m_vertices;
}

template <typename Vector2Container>
const BoundaryTolerance& PolygonBoundaryCheck<Vector2Container>::tolerance()
    const {
  return m_tolerance;
}

template <typename Vector2Container>
bool PolygonBoundaryCheck<Vector2Container>::inside(
    const Vector2& point,
    const std::optional<SquareMatrix2>& jacobianOpt) const {
  if (m_tolerance.isInfinite()) {
    // The null boundary check always succeeds
    return true;
  }

  if (detail::VerticesHelper::isInsidePolygon(point, m_vertices)) {
    // If the point falls inside the polygon, the check always succeeds
    return true;
  }

  if (!m_tolerance.hasTolerance()) {
    // Outside of the polygon, since we've eliminated the case of an absence of
    // check above, we know we'll always fail if the tolerance is zero.
    //
    // This allows us to avoid the expensive computeClosestPointOnPolygon
    // computation in this simple case.
    return false;
  }

  // TODO: When tolerance is not 0, we could also avoid this computation in
  //       some cases by testing against a bounding box of the polygon, padded
  //       on each side with our tolerance. Check if this optimization is
  //       worthwhile in some production workflows, and if so implement it.

  SquareMatrix2 metric = m_tolerance.getMetric(jacobianOpt);

  // We are outside of the polygon, but there is a tolerance. Must find what
  // the closest point on the polygon is and check if it's within tolerance.
  auto closestPoint = detail::VerticesHelper::computeClosestPointOnPolygon(
      point, m_vertices, metric);

  Vector2 distance = closestPoint - point;

  return m_tolerance.isTolerated(distance, jacobianOpt);
}

}  // namespace Acts
