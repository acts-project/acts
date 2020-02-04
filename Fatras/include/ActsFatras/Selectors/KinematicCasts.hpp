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

namespace ActsFatras {
namespace Casts {

/// Retrieve the transverse absolute distance of the vertex to the origin.
struct Vrho {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return Acts::VectorHelpers::perp(particle.position());
  }
};

/// Retrieve the longitudinal distance of the vertex to the origin.
struct Vz {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.position().z();
  }
};

/// Retrieve the longitudinal absolute distance of the vertex to the origin.
struct AbsVz {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return std::abs(particle.position().z());
  }
};

/// Retrieve the direction pseudo-rapidity.
struct Eta {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return Acts::VectorHelpers::eta(particle.direction());
  }
};

/// Retrieve the direction absolute pseudo-rapidity.
struct AbsEta {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return std::abs(Acts::VectorHelpers::eta(particle.direction()));
  }
};

/// Retrieve the transverse momentum.
struct Pt {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.momentum() *
           Acts::VectorHelpers::perp(particle.direction());
  }
};

/// Retrieve the absolute momentum.
struct P {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.momentum();
  }
};

/// Retrieve the total energy.
struct E {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.energy();
  }
};

}  // namespace Casts
}  // namespace ActsFatras
