// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Charge.hpp"
#include "Acts/EventData/GenericParticleHypothesis.hpp"

namespace Acts {

using SinglyChargedParticleHypothesis =
    GenericParticleHypothesis<SinglyCharged>;
using NeutralParticleHypothesis = GenericParticleHypothesis<Neutral>;
using ParticleHypothesis = GenericParticleHypothesis<AnyCharge>;

}  // namespace Acts
