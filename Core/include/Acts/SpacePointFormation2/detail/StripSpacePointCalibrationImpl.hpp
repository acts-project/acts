// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"
#include "Acts/Utilities/detail/StdSpanLinalg.hpp"

#include <cmath>
#include <span>

namespace Acts::detail {

// Intentionally not using Eigen here as there is a noticeable performance
// (10-20%) penalty for small vector sizes.

inline OuterStripSpacePointCalibrationDetailsDerived
deriveOuterStripSpacePointCalibrationDetails(
    const std::span<const float, 3> ihv, const std::span<const float, 3> ohv,
    const std::span<const float, 3> iosv, const std::span<const float, 3> oc) {
  OuterStripSpacePointCalibrationDetailsDerived result{};
  result.innerCrossOuterHalfVector = stdSpanCross(ihv, ohv);
  result.innerToOuterSeparationCrossOuterHalfVector = stdSpanCross(iosv, ohv);
  result.innerToOuterSeparationCrossInnerHalfVector = stdSpanCross(iosv, ihv);
  result.outerCenter = stdSpanCopy(oc);
  result.outerHalfVector = stdSpanCopy(ohv);
  return result;
}

inline OuterStripSpacePointCalibrationDetailsDerived
deriveOuterStripSpacePointCalibrationDetails(
    const OuterStripSpacePointCalibrationDetails& sp) {
  return deriveOuterStripSpacePointCalibrationDetails(
      sp.innerHalfVector, sp.outerHalfVector, sp.innerToOuterSeparation,
      sp.outerCenter);
}

inline bool calibrateOuterStripSpacePoint(
    const std::span<const float, 3> direction,
    const std::span<const float, 3> ihvCrossOhv,
    const std::span<const float, 3> iosvCrossOhv,
    const std::span<const float, 3> iosvCrossIhv,
    const std::span<const float, 3> oc, const std::span<const float, 3> ohv,
    const std::span<float, 3> calibrated, const float tolerance) {
  // scale = innerStripHalfVector dot (outerStripHalfVector cross direction)
  const float scale = stdSpanDot(direction, ihvCrossOhv);

  // sInner = innerToOuterSeparationVector dot (outerStripHalfVector cross
  // direction) Check if direction is inside the inner detector element
  const float sInner = stdSpanDot(direction, iosvCrossOhv);
  if (std::abs(sInner) > std::abs(scale) * tolerance) {
    return false;
  }

  // sOuter = innerToOuterSeparationVector dot (innerStripHalfVector cross
  // direction) Check if direction is inside the outer detector element
  const float sOuter = stdSpanDot(direction, iosvCrossIhv);
  if (std::abs(sOuter) > std::abs(scale) * tolerance) {
    return false;
  }

  // Corrected position using the outer strip center and direction
  const float sOuterNorm = sOuter / scale;
  stdSpanAddScaled(oc, ohv, sOuterNorm, calibrated);
  return true;
}

inline bool calibrateOuterStripSpacePoint(
    const std::span<const float, 3> direction,
    const OuterStripSpacePointCalibrationDetailsDerived& sp,
    const std::span<float, 3> calibrated, const float tolerance) {
  return detail::calibrateOuterStripSpacePoint(
      direction, sp.innerCrossOuterHalfVector,
      sp.innerToOuterSeparationCrossOuterHalfVector,
      sp.innerToOuterSeparationCrossInnerHalfVector, sp.outerCenter,
      sp.outerHalfVector, calibrated, tolerance);
}

}  // namespace Acts::detail
