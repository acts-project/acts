// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <cmath>

namespace ActsFatras::Casts {

/// Retrieve the transverse absolute distance of the position to the origin.
struct Vrho {
  double operator()(const Particle& particle) const {
    auto particleState = particle.lastState();
    return Acts::VectorHelpers::perp(particleState.position());
  }
};

/// Retrieve the longitudinal distance of the position to the origin.
struct Vz {
  double operator()(const Particle& particle) const {
    auto particleState = particle.lastState();
    return particleState.position().z();
  }
};

/// Retrieve the longitudinal absolute distance of the position to the origin.
struct AbsVz {
  double operator()(const Particle& particle) const {
    auto particleState = particle.lastState();
    return std::abs(particleState.position().z());
  }
};

/// Retrieve the direction pseudo-rapidity.
struct Eta {
  double operator()(const Particle& particle) const {
    auto particleState = particle.lastState();
    return particleState.eta();
  }
};

/// Retrieve the direction absolute pseudo-rapidity.
struct AbsEta {
  double operator()(const Particle& particle) const {
    auto particleState = particle.lastState();
    return std::abs(particleState.eta());
  }
};

/// Retrieve the transverse momentum.
struct Pt {
  double operator()(const Particle& particle) const {
    auto particleState = particle.lastState();
    return particleState.transverseMomentum();
  }
};

/// Retrieve the absolute momentum.
struct P {
  double operator()(const Particle& particle) const {
    auto particleState = particle.lastState();
    return particleState.absoluteMomentum();
  }
};

/// Retrieve the total energy.
struct E {
  double operator()(const Particle& particle) const {
    auto particleState = particle.lastState();
    return particleState.energy();
  }
};

}  // namespace ActsFatras::Casts
