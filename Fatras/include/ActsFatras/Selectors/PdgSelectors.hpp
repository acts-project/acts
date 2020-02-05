// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdlib>

#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {

/// Select particles of one specific type.
///
/// Particle and Antiparticle are treated as two separate types.
template <int Pdg>
struct PdgSelector {
  bool operator()(const Particle &particle) const {
    return (static_cast<int>(particle.pdg()) == Pdg);
  }
};

/// Select particles and antiparticles of one specific type.
template <int Pdg>
struct AbsPdgSelector {
  bool operator()(const Particle &particle) const {
    return (std::abs(static_cast<int>(particle.pdg())) == std::abs(Pdg));
  }
};

/// Select all particles except one specific type.
///
/// Particle and Antiparticle are treated as two separate types.
template <int Pdg>
struct PdgExcluder {
  bool operator()(const Particle &particle) const {
    return (static_cast<int>(particle.pdg()) != Pdg);
  }
};

/// Select all particles except for (anti-)particles of one specific type.
template <int Pdg>
struct AbsPdgExcluder {
  bool operator()(const Particle &particle) const {
    return (std::abs(static_cast<int>(particle.pdg())) != std::abs(Pdg));
  }
};

}  // namespace ActsFatras
