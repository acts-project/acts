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

// The stopping-power parameterizations have discontinuities at these branch
// boundaries. Keep their locations shared with the dense transport Jacobian.
inline constexpr float BetheDensityEffectBetaGammaThreshold = 10.f;
inline constexpr float MuonRadiativeLowEnergyThreshold = 8 * UnitConstants::GeV;
inline constexpr float MuonRadiativeHighEnergyThreshold = UnitConstants::TeV;

}  // namespace Acts::detail
