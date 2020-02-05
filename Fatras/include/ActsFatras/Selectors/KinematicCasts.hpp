// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

#include "Acts/Utilities/Helpers.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {
namespace Casts {

/// Retrieve the transverse absolute distance of the vertex to the origin.
struct Vrho {
  double operator()(const Particle& particle) const {
    return Acts::VectorHelpers::perp(particle.position());
  }
};

/// Retrieve the longitudinal distance of the vertex to the origin.
struct Vz {
  double operator()(const Particle& particle) const {
    return particle.position().z();
  }
};

/// Retrieve the longitudinal absolute distance of the vertex to the origin.
struct AbsVz {
  double operator()(const Particle& particle) const {
    return std::abs(particle.position().z());
  }
};

/// Retrieve the direction pseudo-rapidity.
struct Eta {
  double operator()(const Particle& particle) const {
    return Acts::VectorHelpers::eta(particle.direction());
  }
};

/// Retrieve the direction absolute pseudo-rapidity.
struct AbsEta {
  double operator()(const Particle& particle) const {
    return std::abs(Acts::VectorHelpers::eta(particle.direction()));
  }
};

/// Retrieve the transverse momentum.
struct Pt {
  double operator()(const Particle& particle) const {
    return particle.momentum() *
           Acts::VectorHelpers::perp(particle.direction());
  }
};

/// Retrieve the absolute momentum.
struct P {
  double operator()(const Particle& particle) const {
    return particle.momentum();
  }
};

/// Retrieve the total energy.
struct E {
  double operator()(const Particle& particle) const {
    return particle.energy();
  }
};

}  // namespace Casts
}  // namespace ActsFatras
