// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/SeedConfirmationRange.hpp"

// System include(s).
#include <cstddef>

namespace Acts {
struct SeedFilterConfig {
  // the allowed delta between two inverted seed radii for them to be considered
  // compatible.
  float deltaInvHelixDiameter = 0.00003 * 1. / Acts::UnitConstants::mm;
  // the impact parameters (d0) is multiplied by this factor and subtracted from
  // weight
  float impactWeightFactor = 1.;
  // seed weight increased by this value if a compatible seed has been found.
  float compatSeedWeight = 200.;
  // minimum distance between compatible seeds to be considered for weight boost
  float deltaRMin = 5. * Acts::UnitConstants::mm;
  // in dense environments many seeds may be found per middle space point.
  // only seeds with the highest weight will be kept if this limit is reached.
  unsigned int maxSeedsPerSpM = 10;
  // how often do you want to increase the weight of a seed for finding a
  // compatible seed?
  size_t compatSeedLimit = 2;
  // Tool to apply experiment specific cuts on collected middle space points

  // sort vectors vectors by curvature
  bool curvatureSortingInFilter = false;

  // increment in seed weight if number of compatible seeds is larger than
  // numSeedIncrement
  float seedWeightIncrement = 0;
  float numSeedIncrement = std::numeric_limits<float>::infinity();

  // seed confirmation
  bool seedConfirmation = false;
  // contains parameters for central seed confirmation
  SeedConfirmationRange centralSeedConfirmationRange;
  // contains parameters for forward seed confirmation
  SeedConfirmationRange forwardSeedConfirmationRange;
  // minimum radius for bottom SP in seed confirmation
  float seedConfMinBottomRadius = 60. * Acts::UnitConstants::mm;
  // maximum zOrigin in seed confirmation
  float seedConfMaxZOrigin = 150. * Acts::UnitConstants::mm;
  // minimum impact parameter for seed confirmation
  float minImpactSeedConf = 1. * Acts::UnitConstants::mm;

  // maximum number of lower quality seeds in seed confirmation
  float maxSeedsPerSpMConf = std::numeric_limits<float>::infinity();
  // maximum number of quality seeds in seed confirmation
  float maxQualitySeedsPerSpMConf = std::numeric_limits<float>::infinity();

  // use deltaR instead of top radius
  bool useDeltaRorTopRadius = false;

  SeedFilterConfig toInternalUnits() const {
    using namespace Acts::UnitLiterals;
    SeedFilterConfig config = *this;
    config.deltaRMin /= 1_mm;
    config.deltaInvHelixDiameter /= 1. / 1_mm;
    config.seedConfMinBottomRadius /= 1_mm;
    config.seedConfMaxZOrigin /= 1_mm;
    config.minImpactSeedConf /= 1_mm;

    return config;
  }
};

}  // namespace Acts
