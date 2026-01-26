// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/PdgParticle.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {

/// No-op particle selector that selects all particles.
struct EveryParticle {
  /// Select all particles unconditionally
  /// @param particle The particle to evaluate (unused)
  /// @return Always true
  bool operator()(const Particle &particle) const {
    static_cast<void>(particle);
    return true;
  }
};

/// Select neutral particles.
struct NeutralSelector {
  /// Check if particle is neutral
  /// @param particle The particle to evaluate
  /// @return true if particle charge is zero
  bool operator()(const Particle &particle) const {
    return (particle.charge() == 0.);
  }
};

/// Select all charged particles.
struct ChargedSelector {
  /// Check if particle is charged
  /// @param particle The particle to evaluate
  /// @return true if particle charge is non-zero
  bool operator()(const Particle &particle) const {
    return (particle.charge() != 0.);
  }
};

/// Select positively charged particles.
struct PositiveSelector {
  /// Check if particle is positively charged
  /// @param particle The particle to evaluate
  /// @return true if particle charge is positive
  bool operator()(const Particle &particle) const {
    return (0. < particle.charge());
  }
};

/// Select negatively charged particles.
struct NegativeSelector {
  /// Check if particle is negatively charged
  /// @param particle The particle to evaluate
  /// @return true if particle charge is negative
  bool operator()(const Particle &particle) const {
    return (particle.charge() < 0.);
  }
};

/// Select particles of one specific type.
///
/// Particle and Antiparticle are treated as two separate types.
template <Acts::PdgParticle Pdg>
struct PdgSelector {
  /// Check if particle matches the specific PDG type
  /// @param particle The particle to evaluate
  /// @return true if particle PDG matches the template parameter
  bool operator()(const Particle &particle) const {
    return (particle.pdg() == Pdg);
  }
};

/// Select particles and antiparticles of one specific type.
template <Acts::PdgParticle Pdg>
struct AbsPdgSelector {
  /// Check if particle matches the specific PDG type (ignoring sign)
  /// @param particle The particle to evaluate
  /// @return true if absolute PDG type matches the template parameter
  bool operator()(const Particle &particle) const {
    return (makeAbsolutePdgParticle(particle.pdg()) ==
            makeAbsolutePdgParticle(Pdg));
  }
};

/// Select all particles except one specific type.
///
/// Particle and Antiparticle are treated as two separate types.
template <Acts::PdgParticle Pdg>
struct PdgExcluder {
  /// Check if particle does not match the specific PDG type
  /// @param particle The particle to evaluate
  /// @return true if particle PDG does not match the template parameter
  bool operator()(const Particle &particle) const {
    return (particle.pdg() != Pdg);
  }
};

/// Select all particles except for (anti-)particles of one specific type.
template <Acts::PdgParticle Pdg>
struct AbsPdgExcluder {
  /// Check if particle does not match the specific PDG type (ignoring sign)
  /// @param particle The particle to evaluate
  /// @return true if absolute PDG type does not match the template parameter
  bool operator()(const Particle &particle) const {
    return (makeAbsolutePdgParticle(particle.pdg()) !=
            makeAbsolutePdgParticle(Pdg));
  }
};

}  // namespace ActsFatras
