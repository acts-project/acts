// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"

#include <span>

namespace Acts::detail {

inline StripSpacePointCalibrationDetailsDerived
deriveStripSpacePointCalibrationDetails(
    const StripSpacePointCalibrationDetails& sp) {
  const auto& oiv = sp.outerToInnerGapVector;
  const auto& ihv = sp.innerHalfVector;
  const auto& ohv = sp.outerHalfVector;

  return StripSpacePointCalibrationDetailsDerived{
      .outerToInnerGapCrossInnerHalfVector = {oiv[1] * ihv[2] - oiv[2] * ihv[1],
                                              oiv[2] * ihv[0] - oiv[0] * ihv[2],
                                              oiv[0] * ohv[1] -
                                                  oiv[1] * ohv[0]},
      .outerToInnerGapCrossOuterHalfVector = {oiv[1] * ohv[2] - oiv[2] * ohv[1],
                                              oiv[2] * ohv[0] - oiv[0] * ohv[2],
                                              oiv[0] * ohv[1] -
                                                  oiv[1] * ohv[0]},
      .innerCrossOuterHalfVector = {ihv[1] * ohv[2] - ihv[2] * ohv[1],
                                    ihv[2] * ohv[0] - ihv[0] * ohv[2],
                                    ihv[0] * ohv[1] - ihv[1] * ohv[0]},
      .outerCenter = sp.outerCenter,
      .outerHalfVector = ohv,
  };
}

inline bool calibrateStripSpacePoint(
    const StripSpacePointCalibrationDetailsDerived& sp,
    std::span<const float, 3> direction, std::span<float, 3> calibrated,
    float tolerance) {
  const auto& ihvCrossOhv = sp.innerCrossOuterHalfVector;
  const auto& oivCrossOhv = sp.outerToInnerGapCrossOuterHalfVector;
  const auto& oivCrossIhv = sp.outerToInnerGapCrossInnerHalfVector;

  // scale = innerStripHalfVector dot (outerStripHalfVector cross direction)
  const float scale = direction[0] * ihvCrossOhv[0] +
                      direction[1] * ihvCrossOhv[1] +
                      direction[2] * ihvCrossOhv[2];

  // sInner = outerToInnerGapVector dot (outerStripHalfVector cross direction)
  // Check if direction is inside the inner detector element
  const float sInner = direction[0] * oivCrossOhv[0] +
                       direction[1] * oivCrossOhv[1] +
                       direction[2] * oivCrossOhv[2];
  if (std::abs(sInner) > std::abs(scale) * tolerance) {
    return false;
  }

  // sOuter = outerToInnerGapVector dot (innerStripHalfVector cross direction)
  // Check if direction is inside the outer detector element
  const float sOuter = direction[0] * oivCrossIhv[0] +
                       direction[1] * oivCrossIhv[1] +
                       direction[2] * oivCrossIhv[2];
  if (std::abs(sOuter) > std::abs(scale) * tolerance) {
    return false;
  }

  const auto& oc = sp.outerCenter;
  const auto& ohv = sp.outerHalfVector;

  // Corrected position using the outer strip center and direction
  const float sOuterNorm = sOuter / scale;
  calibrated[0] = oc[0] + ohv[0] * sOuterNorm;
  calibrated[1] = oc[1] + ohv[1] * sOuterNorm;
  calibrated[2] = oc[2] + ohv[2] * sOuterNorm;
  return true;
}

}  // namespace Acts::detail
