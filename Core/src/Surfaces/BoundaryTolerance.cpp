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

namespace Acts {

BoundaryTolerance::ToleranceMode BoundaryTolerance::toleranceMode() const {
  using enum ToleranceMode;

  if (isInfinite()) {
    return Extend;
  }

  if (isNone()) {
    return None;
  }

  if (const auto* absoluteEuclidean = getVariantPtr<AbsoluteEuclideanParams>();
      absoluteEuclidean != nullptr) {
    if (absoluteEuclidean->tolerance == 0) {
      return None;
    } else if (absoluteEuclidean->tolerance > 0) {
      return Extend;
    } else {
      return Shrink;
    }
  }

  if (const auto* chi2Bound = getVariantPtr<Chi2BoundParams>();
      chi2Bound != nullptr) {
    if (chi2Bound->maxChi2 == 0) {
      return None;
    } else if (chi2Bound->maxChi2 >= 0) {
      return Extend;
    } else {
      return Shrink;
    }
  }

  if (const auto* chi2Cartesian = getVariantPtr<Chi2CartesianParams>();
      chi2Cartesian != nullptr) {
    if (chi2Cartesian->maxChi2 == 0) {
      return None;
    } else if (chi2Cartesian->maxChi2 >= 0) {
      return Extend;
    } else {
      return Shrink;
    }
  }

  throw std::logic_error("Unsupported tolerance type");
}

bool BoundaryTolerance::isTolerated(
    const Vector2& boundDelta, const SquareMatrix2& boundToCartesian) const {
  if (isInfinite()) {
    return true;
  }

  if (isNone()) {
    return boundDelta == Vector2::Zero();
  }

  if (const auto* chi2Bound = getVariantPtr<Chi2BoundParams>();
      chi2Bound != nullptr) {
    double chi2 =
        boundDelta.transpose() * chi2Bound->weightMatrix() * boundDelta;
    return std::copysign(chi2, chi2Bound->maxChi2) <= chi2Bound->maxChi2;
  }

  Vector2 cartesianDelta = boundToCartesian * boundDelta;

  if (const auto* absoluteEuclidean = getVariantPtr<AbsoluteEuclideanParams>();
      absoluteEuclidean != nullptr) {
    return std::copysign(cartesianDelta.norm(), absoluteEuclidean->tolerance) <=
           absoluteEuclidean->tolerance;
  }

  if (const auto* chi2Cartesian = getVariantPtr<Chi2CartesianParams>();
      chi2Cartesian != nullptr) {
    double chi2 = cartesianDelta.transpose() * chi2Cartesian->weightMatrix() *
                  cartesianDelta;
    return std::copysign(chi2, chi2Cartesian->maxChi2) <=
           chi2Cartesian->maxChi2;
  }

  throw std::logic_error("Unsupported tolerance type");
}

void BoundaryTolerance::print(std::ostream& ostr) const {
  if (isInfinite()) {
    ostr << "no boundary limit applied";
  } else if (isNone()) {
    ostr << "strict boundary limit";
  } else if (hasAbsoluteEuclidean()) {
    ostr << "boundary limit  dX < " << asAbsoluteEuclidean().tolerance;
  } else if (hasChi2Bound()) {
    ostr << "local chi2 " << asChi2Bound().maxChi2 << " with weights:\n"
         << asChi2Bound().weightMatrix() << "\n";
  } else if (hasChi2Cartesian()) {
    ostr << "local chi2 " << asChi2Cartesian().maxChi2 << " with weights:\n"
         << asChi2Cartesian().weightMatrix() << "\n";
  }
}

}  // namespace Acts
