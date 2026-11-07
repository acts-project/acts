// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Synthetic/EventConfig.hpp"

#include "Acts/Definitions/ParticleData.hpp"

#include <optional>
#include <stdexcept>

namespace ActsFatras::Synthetic {

float EventConfig::particleMass() const {
  const std::optional<float> mass =
      Acts::findMass(Acts::makeAbsolutePdgParticle(particlePdg));
  if (!mass.has_value()) {
    throw std::invalid_argument(
        "EventConfig: the core library has no mass for this particlePdg");
  }
  return *mass;
}

}  // namespace ActsFatras::Synthetic
