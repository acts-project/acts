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
  /// @brief Extract transverse distance from origin
  /// @param particle The particle to extract from
  /// @return Transverse distance rho = sqrt(x² + y²) from particle position
  double operator()(const Particle& particle) const {
    return std::hypot(particle.position().x(), particle.position().y());
  }
};

/// Retrieve the longitudinal distance of the position to the origin.
struct Vz {
  /// @brief Extract longitudinal distance from origin
  /// @param particle The particle to extract from
  /// @return Z-coordinate of particle position
  double operator()(const Particle& particle) const {
    return particle.position().z();
  }
};

/// Retrieve the longitudinal absolute distance of the position to the origin.
struct AbsVz {
  /// @brief Extract absolute longitudinal distance from origin
  /// @param particle The particle to extract from
  /// @return Absolute value of Z-coordinate of particle position
  double operator()(const Particle& particle) const {
    return std::abs(particle.position().z());
  }
};

/// Retrieve the direction pseudo-rapidity.
struct Eta {
  /// @brief Extract direction pseudo-rapidity
  /// @param particle The particle to extract from
  /// @return Pseudo-rapidity η = atanh(p_z/p) from particle direction
  double operator()(const Particle& particle) const {
    // particle direction is always normalized, i.e. dz = pz / p
    return std::atanh(particle.direction().z());
  }
};

/// Retrieve the direction absolute pseudo-rapidity.
struct AbsEta {
  /// @brief Extract absolute direction pseudo-rapidity
  /// @param particle The particle to extract from
  /// @return Absolute pseudo-rapidity |η| = atanh(|p_z|/p) from particle direction
  double operator()(const Particle& particle) const {
    // particle direction is always normalized, i.e. dz = pz / p
    return std::atanh(std::abs(particle.direction().z()));
  }
};

/// Retrieve the transverse momentum.
struct Pt {
  /// @brief Extract transverse momentum
  /// @param particle The particle to extract from
  /// @return Transverse momentum p_T = p × sin(θ) from particle momentum
  double operator()(const Particle& particle) const {
    // particle direction is always normalized, i.e. dt²+dz²=1 w/ dt²=dx²+dy²
    return particle.absoluteMomentum() *
           Acts::VectorHelpers::perp(particle.direction());
  }
};

/// Retrieve the absolute momentum.
struct P {
  /// @brief Extract absolute momentum magnitude
  /// @param particle The particle to extract from
  /// @return Total momentum magnitude |p| of the particle
  double operator()(const Particle& particle) const {
    return particle.absoluteMomentum();
  }
};

/// Retrieve the total energy.
struct E {
  /// @brief Extract total energy
  /// @param particle The particle to extract from
  /// @return Total energy E = sqrt(p²c² + m²c⁴) of the particle
  double operator()(const Particle& particle) const {
    return particle.energy();
  }
};

}  // namespace ActsFatras::Casts
