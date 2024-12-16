// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Surfaces/BoundaryTolerance.hpp"

#include "Acts/Definitions/Algebra.hpp"

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

BoundaryTolerance::ToleranceMode BoundaryTolerance::toleranceMode() const {
  using enum ToleranceMode;
  if (isInfinite()) {
    return Extend;
  }

  if (isNone()) {
    return None;
  }

  if (const auto* absoluteBound = getVariantPtr<AbsoluteBound>();
      absoluteBound != nullptr) {
    if (absoluteBound->tolerance0 == 0. && absoluteBound->tolerance1 == 0.) {
      return None;
    }

    // std::cout << absoluteBound->tolerance0 << " "
    //           << std::copysign(1., absoluteBound->tolerance0) << std::endl;
    // std::cout << absoluteBound->tolerance1 << " "
    //           << std::copysign(1., absoluteBound->tolerance1) << std::endl;

    if (std::copysign(1., absoluteBound->tolerance0) !=
        std::copysign(1., absoluteBound->tolerance1)) {
      throw std::logic_error("Inconsistent tolerance signs are not supported");
    }

    if (absoluteBound->tolerance0 > 0. || absoluteBound->tolerance1 > 0.) {
      return Extend;
    } else {
      return Shrink;
    }
  }

  if (const auto* absoluteCartesian = getVariantPtr<AbsoluteCartesian>();
      absoluteCartesian != nullptr) {
    if (absoluteCartesian->tolerance0 == 0. &&
        absoluteCartesian->tolerance1 == 0.) {
      return None;
    }

    if (std::copysign(1., absoluteCartesian->tolerance0) !=
        std::copysign(1., absoluteCartesian->tolerance1)) {
      throw std::logic_error("Inconsistent tolerance signs are not supported");
    }

    if (absoluteCartesian->tolerance0 > 0. ||
        absoluteCartesian->tolerance1 > 0.) {
      return Extend;
    } else {
      return Shrink;
    }
  }

  if (const auto* absoluteEuclidean = getVariantPtr<AbsoluteEuclidean>();
      absoluteEuclidean != nullptr) {
    if (absoluteEuclidean->tolerance == 0.) {
      return None;
    } else if (absoluteEuclidean->tolerance > 0.) {
      return Extend;
    } else {
      return Shrink;
    }
  }

  if (const auto* chi2Bound = getVariantPtr<Chi2Bound>();
      chi2Bound != nullptr) {
    if (chi2Bound->maxChi2 == 0.) {
      return None;
    } else if (chi2Bound->maxChi2 > 0.) {
      return Extend;
    } else {
      return Shrink;
    }
  }

  assert(false && "Unsupported tolerance type");
  return None;
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

  if (const auto* absoluteBound = getVariantPtr<AbsoluteBound>();
      absoluteBound != nullptr) {
    bool tol0 = false;
    if (absoluteBound->tolerance0 < 0) {
      tol0 = std::abs(distance[0]) > std::abs(absoluteBound->tolerance0);
    } else {
      tol0 = std::abs(distance[0]) <= absoluteBound->tolerance0;
    }

    bool tol1 = false;
    if (absoluteBound->tolerance1 < 0) {
      tol1 = std::abs(distance[1]) > std::abs(absoluteBound->tolerance1);
    } else {
      tol1 = std::abs(distance[1]) <= absoluteBound->tolerance1;
    }

    return tol0 && tol1;
  }

  if (const auto* chi2Bound = getVariantPtr<Chi2Bound>();
      chi2Bound != nullptr) {
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

  if (const auto* absoluteCartesian = getVariantPtr<AbsoluteCartesian>();
      absoluteCartesian != nullptr) {
    bool tol0 = false;

    if (absoluteCartesian->tolerance0 < 0) {
      tol0 = cartesianDistance[0] > std::abs(absoluteCartesian->tolerance0);
    } else {
      tol0 = std::abs(cartesianDistance[0]) <= absoluteCartesian->tolerance0;
    }

    bool tol1 = false;
    if (absoluteCartesian->tolerance1 < 0) {
      tol1 = cartesianDistance[1] > std::abs(absoluteCartesian->tolerance1);
    } else {
      tol1 = std::abs(cartesianDistance[1]) <= absoluteCartesian->tolerance1;
    }

    return tol0 && tol1;
  }

  if (const auto* absoluteEuclidean = getVariantPtr<AbsoluteEuclidean>();
      absoluteEuclidean != nullptr) {
    if (absoluteEuclidean->tolerance < 0) {
      return cartesianDistance.norm() > std::abs(absoluteEuclidean->tolerance);
    } else {
      return cartesianDistance.norm() <= absoluteEuclidean->tolerance;
    }
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

  if (const auto* chi2Bound = getVariantPtr<BoundaryTolerance::Chi2Bound>();
      chi2Bound != nullptr) {
    metric = chi2Bound->weight;
  } else if (!isCartesian) {
    const auto& jacobian = *jacobianOpt;
    metric = jacobian.transpose() * jacobian;
  }

  return metric;
}

}  // namespace Acts
