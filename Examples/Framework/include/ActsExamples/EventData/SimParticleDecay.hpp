// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/SimParticle.hpp"

#include <unordered_map>
#include <vector>

namespace ActsExamples {

/// Generator-selected decay of a transported particle, in ACTS units.
/// Daughter momenta are expressed in the laboratory frame at generation.
struct SimParticleDecay {
  double properTime = 0;
  Acts::Vector4 endPosition = Acts::Vector4::Zero();
  std::vector<SimParticle> daughters;
};

/// Keyed by the generator-assigned ACTS barcode of each decaying particle.
using SimParticleDecayGraph = std::unordered_map<SimBarcode, SimParticleDecay>;

}  // namespace ActsExamples
