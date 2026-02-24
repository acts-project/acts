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
#include "Acts/Utilities/Diagnostics.hpp"

namespace Acts {

ACTS_PUSH_IGNORE_DEPRECATED()
/// @brief Type alias for bound track parameters of singly charged particles
/// @details Track parameters defined in a local coordinate system for particles with charge ±1
using SinglyChargedBoundTrackParameters [[deprecated(
    "Use BoundTrackParameters with one charge magnitude instead")]] =
    GenericBoundTrackParameters<SinglyChargedParticleHypothesis>;
/// Type alias for free track parameters of singly charged particles
using SinglyChargedFreeTrackParameters [[deprecated(
    "Use FreeTrackParameters with one charge magnitude instead")]] =
    GenericFreeTrackParameters<SinglyChargedParticleHypothesis>;

/// @brief Type alias for bound track parameters of neutral particles
/// @details Track parameters for neutral particles in a bound coordinate system
using NeutralBoundTrackParameters [[deprecated(
    "Use BoundTrackParameters with zero charge magnitude instead")]] =
    GenericBoundTrackParameters<NeutralParticleHypothesis>;
/// Type alias for free track parameters of neutral particles
using NeutralFreeTrackParameters [[deprecated(
    "Use FreeTrackParameters with zero charge magnitude instead")]] =
    GenericFreeTrackParameters<NeutralParticleHypothesis>;
ACTS_POP_IGNORE_DEPRECATED()

/// @brief BoundTrackParameters can hold any kind of charge
using BoundTrackParameters = GenericBoundTrackParameters<ParticleHypothesis>;
/// @brief FreeTrackParameters can hold any kind of charge
using FreeTrackParameters = GenericFreeTrackParameters<ParticleHypothesis>;

}  // namespace Acts
