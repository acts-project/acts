// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>

namespace Acts {
/// @brief contains parameters for  seed confirmation
struct SeedConfirmationRangeConfig {
  // z minimum and maximum of middle component of the seed used to define the
  // region of the detector for seed confirmation
  float zMinSeedConf =
      std::numeric_limits<float>::lowest();  // Acts::UnitConstants::mm
  float zMaxSeedConf =
      std::numeric_limits<float>::max();  // Acts::UnitConstants::mm
  // radius of bottom component of seed that is used to define the number of
  // compatible top required
  float rMaxSeedConf =
      std::numeric_limits<float>::max();  // Acts::UnitConstants::mm

  // number of compatible top SPs of seed if bottom radius is larger than
  // rMaxSeedConf
  size_t nTopForLargeR = 0;
  // number of compatible top SPs of seed if bottom radius is smaller than
  // rMaxSeedConf
  size_t nTopForSmallR = 0;

  // minimum radius for bottom SP in seed confirmation
  float seedConfMinBottomRadius = 60. * Acts::UnitConstants::mm;
  // maximum zOrigin in seed confirmation
  float seedConfMaxZOrigin = 150. * Acts::UnitConstants::mm;
  // minimum impact parameter for seed confirmation
  float minImpactSeedConf = 1. * Acts::UnitConstants::mm;
};
}  // namespace Acts
