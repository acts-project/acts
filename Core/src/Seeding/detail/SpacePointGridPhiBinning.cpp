// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/detail/SpacePointGridPhiBinning.hpp"

#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace Acts::detail {

int computeSpacePointGridPhiBins(const SpacePointGridPhiBinningConfig& config,
                                 const Logger& logger) {
  if (config.bFieldInZ == 0) {
    // for no magnetic field, use the maximum number of phi bins
    ACTS_VERBOSE(
        "B-Field is 0 (z-coordinate), setting the number of bins in phi to "
        << config.maxPhiBins);
    return config.maxPhiBins;
  }

  // calculate circle intersections of helix and max detector radius in mm.
  // bFieldInZ is in (pT/radius) natively, no need for conversion
  const float minHelixRadius = config.minPt / config.bFieldInZ;

  // sanity check: if yOuter takes the square root of a negative number
  if (minHelixRadius < config.rMax * 0.5) {
    throw std::domain_error(
        "The value of minHelixRadius cannot be smaller than rMax / 2. Please "
        "check the configuration of bFieldInZ and minPt");
  }

  // x = rMax^2 / (2 * minHelixRadius)
  // y = cathetus(rMax, x)
  // R / x = 2 * minHelixRadius / R
  // outerAngle = x / y = x / sqrt(rMax^2 - x^2) = 1 / cath(rMax / x, 1)
  //            = 1 / cath(2 * minHelixRadius / rMax, 1)
  const float outerAngle =
      std::atan(1.f / fastCathetus(2 * minHelixRadius / config.rMax, 1));
  // intersection of helix and max detector radius minus maximum R distance
  // from middle SP to top SP
  float innerAngle = 0;
  float rMin = config.rMax;
  if (config.rMax > config.deltaRMax) {
    const float innerCircleR = config.rMax - config.deltaRMax;
    rMin = innerCircleR;
    innerAngle =
        std::atan(1.f / fastCathetus(2 * minHelixRadius / innerCircleR, 1));
  }

  // evaluating the azimutal deflection including the maximum impact
  // parameter. A track with |d0| >= r cannot reach radius r, so the azimuthal
  // deflection saturates at pi/2: clamp the asin arguments to [0, 1] so that
  // impactMax >= (rMax - deltaRMax) falls back to coarse / full-2pi phi
  // coverage instead of producing a NaN (and hence a zero phi-bin count and a
  // "Invalid binning" exception).
  const float sinInner = std::min(1.f, config.impactMax / rMin);
  const float sinOuter = std::min(1.f, config.impactMax / config.rMax);
  const float deltaAngleWithMaxD0 =
      std::abs(std::asin(sinInner) - std::asin(sinOuter));

  // evaluating delta Phi based on the inner and outer angle, and the azimutal
  // deflection including the maximum impact parameter
  // Divide by config.phiBinDeflectionCoverage since we combine
  // config.phiBinDeflectionCoverage number of consecutive phi bins in the
  // seed making step. So each individual bin should cover
  // 1/config.phiBinDeflectionCoverage of the maximum expected azimutal
  // deflection
  const float deltaPhi = (outerAngle - innerAngle + deltaAngleWithMaxD0) /
                         config.phiBinDeflectionCoverage;

  // sanity check: if the delta phi is equal to or less than zero, we'll be
  // creating an infinite or a negative number of bins, which would be bad!
  if (deltaPhi <= 0.f) {
    throw std::domain_error(
        "Delta phi value is equal to or less than zero, leading to an "
        "impossible number of bins (negative or infinite)");
  }

  // divide 2pi by angle delta to get number of phi-bins
  // size is always 2pi even for regions of interest
  const int phiBins =
      static_cast<int>(std::ceil(2 * std::numbers::pi / deltaPhi));

  // set protection for large number of bins, by default it is large
  return std::min(phiBins, config.maxPhiBins);
}

}  // namespace Acts::detail
