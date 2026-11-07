// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

namespace Acts::detail {

// The stopping power jumps at these thresholds. The dense stepper jacobian
// needs the same values.
inline constexpr float BetheDensityEffectBetaGammaThreshold = 10.f;
// TODO magic number 8 GeV
inline constexpr float MuonRadiativeLowEnergyThreshold = 8 * UnitConstants::GeV;
inline constexpr float MuonRadiativeHighEnergyThreshold = UnitConstants::TeV;

}  // namespace Acts::detail
