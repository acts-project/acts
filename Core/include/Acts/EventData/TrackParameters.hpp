// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/GenericCurvilinearTrackParameters.hpp"
#include "Acts/EventData/GenericFreeTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"

namespace Acts {

using SinglyChargedBoundTrackParameters =
    GenericBoundTrackParameters<SinglyChargedParticleHypothesis>;
using SinglyChargedCurvilinearTrackParameters =
    GenericCurvilinearTrackParameters<SinglyChargedParticleHypothesis>;
using SinglyChargedFreeTrackParameters =
    GenericFreeTrackParameters<SinglyChargedParticleHypothesis>;

using NeutralBoundTrackParameters =
    GenericBoundTrackParameters<NeutralParticleHypothesis>;
using NeutralCurvilinearTrackParameters =
    GenericCurvilinearTrackParameters<NeutralParticleHypothesis>;
using NeutralFreeTrackParameters =
    GenericFreeTrackParameters<NeutralParticleHypothesis>;

using BoundTrackParameters = GenericBoundTrackParameters<ParticleHypothesis>;
using CurvilinearTrackParameters =
    GenericCurvilinearTrackParameters<ParticleHypothesis>;
using FreeTrackParameters = GenericFreeTrackParameters<ParticleHypothesis>;

}  // namespace Acts
