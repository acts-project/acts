// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/BoundaryCheck.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <iostream>
#include <stdexcept>
#include <utility>

namespace Acts {

BoundaryTolerance::BoundaryTolerance(const Infinite& infinite)
    : m_variant{infinite} {}

BoundaryTolerance::BoundaryTolerance(const None& none) : m_variant{none} {}

BoundaryTolerance::BoundaryTolerance(const AbsoluteBound& absoluteBound)
    : m_variant{absoluteBound} {}

BoundaryTolerance::BoundaryTolerance(const AbsoluteCartesian& absoluteCartesian)
    : m_variant{absoluteCartesian} {}

BoundaryTolerance::BoundaryTolerance(const AbsoluteEuclidean& absoluteEuclidean)
    : m_variant{absoluteEuclidean} {}

BoundaryTolerance::BoundaryTolerance(const Chi2Bound& chi2Bound)
    : m_variant{chi2Bound} {}

BoundaryTolerance::BoundaryTolerance(Variant variant)
    : m_variant{std::move(variant)} {}

bool BoundaryTolerance::isInfinite() const {
  return holdsVariant<Infinite>();
}

bool BoundaryTolerance::isNone() const {
  return holdsVariant<None>();
}

bool BoundaryTolerance::hasAbsoluteBound(bool isCartesian) const {
  return holdsVariant<None>() || holdsVariant<AbsoluteBound>() ||
         (isCartesian && holdsVariant<AbsoluteCartesian>());
}

bool BoundaryTolerance::hasAbsoluteCartesian() const {
  return holdsVariant<AbsoluteCartesian>();
}

bool BoundaryTolerance::hasAbsoluteEuclidean() const {
  return holdsVariant<AbsoluteEuclidean>();
}

bool BoundaryTolerance::hasChi2Bound() const {
  return holdsVariant<Chi2Bound>();
}

bool BoundaryTolerance::hasTolerance() const {
  if (isInfinite()) {
    return true;
  }

  if (isNone()) {
    return false;
  }

  if (auto absoluteBound = getVariantPtr<AbsoluteBound>();
      absoluteBound != nullptr) {
    return absoluteBound->tolerance0 != 0. || absoluteBound->tolerance1 != 0.;
  }

  if (auto absoluteCartesian = getVariantPtr<AbsoluteCartesian>();
      absoluteCartesian != nullptr) {
    return absoluteCartesian->tolerance0 != 0. ||
           absoluteCartesian->tolerance1 != 0.;
  }

  if (auto absoluteEuclidean = getVariantPtr<AbsoluteEuclidean>();
      absoluteEuclidean != nullptr) {
    return absoluteEuclidean->tolerance != 0.;
  }

  if (auto chi2Bound = getVariantPtr<Chi2Bound>(); chi2Bound != nullptr) {
    return chi2Bound->maxChi2 != 0.;
  }

  assert(false && "Unsupported tolerance type");
  return false;
}

BoundaryTolerance::AbsoluteBound BoundaryTolerance::asAbsoluteBound(
    bool isCartesian) const {
  if (isNone()) {
    return AbsoluteBound{0., 0.};
  }

  if (isCartesian && hasAbsoluteCartesian()) {
    const auto& cartesian = getVariant<AbsoluteCartesian>();
    return AbsoluteBound{cartesian.tolerance0, cartesian.tolerance1};
  }

  return getVariant<AbsoluteBound>();
}

const BoundaryTolerance::AbsoluteCartesian&
BoundaryTolerance::asAbsoluteCartesian() const {
  return getVariant<AbsoluteCartesian>();
}

const BoundaryTolerance::AbsoluteEuclidean&
BoundaryTolerance::asAbsoluteEuclidean() const {
  return getVariant<AbsoluteEuclidean>();
}

const BoundaryTolerance::Chi2Bound& BoundaryTolerance::asChi2Bound() const {
  return getVariant<Chi2Bound>();
}

std::optional<BoundaryTolerance::AbsoluteBound>
BoundaryTolerance::asAbsoluteBoundOpt(bool isCartesian) const {
  return hasAbsoluteBound(isCartesian)
             ? std::optional(asAbsoluteBound(isCartesian))
             : std::nullopt;
}

bool BoundaryTolerance::isTolerated(
    const Vector2& distance,
    const std::optional<SquareMatrix2>& jacobianOpt) const {
  if (isInfinite()) {
    return true;
  }

  if (isNone()) {
    return distance == Vector2::Zero();
  }

  if (auto absoluteBound = getVariantPtr<AbsoluteBound>();
      absoluteBound != nullptr) {
    return std::abs(distance[0]) <= absoluteBound->tolerance0 &&
           std::abs(distance[1]) <= absoluteBound->tolerance1;
  }

  if (auto chi2Bound = getVariantPtr<Chi2Bound>(); chi2Bound != nullptr) {
    double chi2 = distance.transpose() * chi2Bound->weight * distance;
    // Mahalanobis distances mean is 2 in 2-dim. cut is 1-d sigma.
    return chi2 <= 2 * chi2Bound->maxChi2;
  }

  bool isCartesian = !jacobianOpt.has_value();
  Vector2 cartesianDistance;
  if (isCartesian) {
    cartesianDistance = distance;
  } else {
    const auto& jacobian = *jacobianOpt;
    cartesianDistance = jacobian * distance;
  }

  if (auto absoluteCartesian = getVariantPtr<AbsoluteCartesian>();
      absoluteCartesian != nullptr) {
    return std::abs(cartesianDistance[0]) <= absoluteCartesian->tolerance0 &&
           std::abs(cartesianDistance[1]) <= absoluteCartesian->tolerance1;
  }

  if (auto absoluteEuclidean = getVariantPtr<AbsoluteEuclidean>();
      absoluteEuclidean != nullptr) {
    return cartesianDistance.norm() <= absoluteEuclidean->tolerance;
  }

  throw std::logic_error("Unsupported tolerance type");
}

bool BoundaryTolerance::hasMetric(bool hasJacobian) const {
  return hasJacobian || hasChi2Bound();
}

SquareMatrix2 BoundaryTolerance::getMetric(
    const std::optional<SquareMatrix2>& jacobianOpt) const {
  bool isCartesian = !jacobianOpt.has_value();
  SquareMatrix2 metric = SquareMatrix2::Identity();

  if (auto chi2Bound = getVariantPtr<BoundaryTolerance::Chi2Bound>();
      chi2Bound != nullptr) {
    metric = chi2Bound->weight;
  } else if (!isCartesian) {
    const auto& jacobian = *jacobianOpt;
    metric = jacobian.transpose() * jacobian;
  }

  return metric;
}

AlignedBoxBoundaryCheck::AlignedBoxBoundaryCheck(const Vector2& lowerLeft,
                                                 const Vector2& upperRight,
                                                 BoundaryTolerance tolerance)
    : m_lowerLeft(lowerLeft),
      m_upperRight(upperRight),
      m_tolerance(std::move(tolerance)) {}

const Vector2& AlignedBoxBoundaryCheck::lowerLeft() const {
  return m_lowerLeft;
}

const Vector2& AlignedBoxBoundaryCheck::upperRight() const {
  return m_upperRight;
}

std::array<Vector2, 4> AlignedBoxBoundaryCheck::vertices() const {
  return {{m_lowerLeft,
           {m_upperRight[0], m_lowerLeft[1]},
           m_upperRight,
           {m_lowerLeft[0], m_upperRight[1]}}};
}

const BoundaryTolerance& AlignedBoxBoundaryCheck::tolerance() const {
  return m_tolerance;
}

bool AlignedBoxBoundaryCheck::inside(
    const Vector2& point,
    const std::optional<SquareMatrix2>& jacobianOpt) const {
  if (m_tolerance.isInfinite()) {
    return true;
  }

  if (detail::VerticesHelper::isInsideRectangle(point, m_lowerLeft,
                                                m_upperRight)) {
    return true;
  }

  if (!m_tolerance.hasTolerance()) {
    return false;
  }

  Vector2 closestPoint;

  if (!m_tolerance.hasMetric(jacobianOpt.has_value())) {
    closestPoint =
        detail::VerticesHelper::computeEuclideanClosestPointOnRectangle(
            point, m_lowerLeft, m_upperRight);
  } else {
    // TODO there might be a more optimal way to compute the closest point to a
    // box with metric

    SquareMatrix2 metric = m_tolerance.getMetric(jacobianOpt);

    closestPoint = detail::VerticesHelper::computeClosestPointOnPolygon(
        point, vertices(), metric);
  }

  Vector2 distance = closestPoint - point;

  return m_tolerance.isTolerated(distance, jacobianOpt);
}

}  // namespace Acts
