// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "ACTFW/EventData/IndexContainers.hpp"
#include "ACTFW/EventData/ProtoTrack.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

namespace FW {

/// Associate a particle to its hit count within a proto track.
struct ParticleHitCount {
  ActsFatras::Barcode particleId;
  std::size_t hitCount;
};

/// Identify all particles that contribute to the proto track and count hits.
void identifyContributingParticles(
    const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
    const ProtoTrack& protoTrack,
    std::vector<ParticleHitCount>& particleHitCount);

}  // namespace FW
