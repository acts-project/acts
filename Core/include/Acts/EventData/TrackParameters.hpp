// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/GenericFreeTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"

namespace Acts {

/// @brief Type alias for bound track parameters of singly charged particles
/// @details Track parameters defined in a local coordinate system for particles with charge ±1
using SinglyChargedBoundTrackParameters =
    GenericBoundTrackParameters<SinglyChargedParticleHypothesis>;
/// Type alias for free track parameters of singly charged particles
using SinglyChargedFreeTrackParameters =
    GenericFreeTrackParameters<SinglyChargedParticleHypothesis>;

/// @brief Type alias for bound track parameters of neutral particles
/// @details Track parameters for neutral particles in a bound coordinate system
using NeutralBoundTrackParameters =
    GenericBoundTrackParameters<NeutralParticleHypothesis>;
/// Type alias for free track parameters of neutral particles
using NeutralFreeTrackParameters =
    GenericFreeTrackParameters<NeutralParticleHypothesis>;

/// @brief BoundTrackParameters can hold any kind of charge
using BoundTrackParameters = GenericBoundTrackParameters<ParticleHypothesis>;
/// @brief FreeTrackParameters can hold any kind of charge
using FreeTrackParameters = GenericFreeTrackParameters<ParticleHypothesis>;

}  // namespace Acts
