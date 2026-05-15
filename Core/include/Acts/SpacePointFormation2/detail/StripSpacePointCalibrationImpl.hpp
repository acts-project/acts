// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"

#include <cmath>

namespace Acts::detail {

inline OuterStripSpacePointCalibrationDetailsDerived
deriveOuterStripSpacePointCalibrationDetails(const std::array<float, 3> ihv,
                                             const std::array<float, 3> ohv,
                                             const std::array<float, 3> iosv,
                                             const std::array<float, 3> oc) {
  OuterStripSpacePointCalibrationDetailsDerived result{};
  result.innerCrossOuterHalfVector = {ihv[1] * ohv[2] - ihv[2] * ohv[1],
                                      ihv[2] * ohv[0] - ihv[0] * ohv[2],
                                      ihv[0] * ohv[1] - ihv[1] * ohv[0]};
  result.innerToOuterSeparationCrossOuterHalfVector = {
      iosv[1] * ohv[2] - iosv[2] * ohv[1], iosv[2] * ohv[0] - iosv[0] * ohv[2],
      iosv[0] * ohv[1] - iosv[1] * ohv[0]};
  result.innerToOuterSeparationCrossInnerHalfVector = {
      iosv[1] * ihv[2] - iosv[2] * ihv[1], iosv[2] * ihv[0] - iosv[0] * ihv[2],
      iosv[0] * ihv[1] - iosv[1] * ihv[0]};
  result.outerCenter = oc;
  result.outerHalfVector = ohv;
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
    const std::array<float, 3>& direction,
    const std::array<float, 3>& ihvCrossOhv,
    const std::array<float, 3>& iosvCrossOhv,
    const std::array<float, 3>& iosvCrossIhv, const std::array<float, 3>& oc,
    const std::array<float, 3>& ohv, std::array<float, 3>& calibrated,
    float tolerance) {
  // scale = innerStripHalfVector dot (outerStripHalfVector cross direction)
  const float scale = direction[0] * ihvCrossOhv[0] +
                      direction[1] * ihvCrossOhv[1] +
                      direction[2] * ihvCrossOhv[2];

  // sInner = innerToOuterSeparationVector dot (outerStripHalfVector cross
  // direction) Check if direction is inside the inner detector element
  const float sInner = direction[0] * iosvCrossOhv[0] +
                       direction[1] * iosvCrossOhv[1] +
                       direction[2] * iosvCrossOhv[2];
  if (std::abs(sInner) > std::abs(scale) * tolerance) {
    return false;
  }

  // sOuter = innerToOuterSeparationVector dot (innerStripHalfVector cross
  // direction) Check if direction is inside the outer detector element
  const float sOuter = direction[0] * iosvCrossIhv[0] +
                       direction[1] * iosvCrossIhv[1] +
                       direction[2] * iosvCrossIhv[2];
  if (std::abs(sOuter) > std::abs(scale) * tolerance) {
    return false;
  }

  // Corrected position using the outer strip center and direction
  const float sOuterNorm = sOuter / scale;
  calibrated[0] = oc[0] + ohv[0] * sOuterNorm;
  calibrated[1] = oc[1] + ohv[1] * sOuterNorm;
  calibrated[2] = oc[2] + ohv[2] * sOuterNorm;
  return true;
}

inline bool calibrateOuterStripSpacePoint(
    const std::array<float, 3>& direction,
    const OuterStripSpacePointCalibrationDetailsDerived& sp,
    std::array<float, 3>& calibrated, const float tolerance) {
  return detail::calibrateOuterStripSpacePoint(
      direction, sp.innerCrossOuterHalfVector,
      sp.innerToOuterSeparationCrossOuterHalfVector,
      sp.innerToOuterSeparationCrossInnerHalfVector, sp.outerCenter,
      sp.outerHalfVector, calibrated, tolerance);
}

}  // namespace Acts::detail
