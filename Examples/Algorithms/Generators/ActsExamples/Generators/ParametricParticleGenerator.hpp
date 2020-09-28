// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/PdgParticle.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace ActsExamples {

/// Generate particles from uniform parameter distributions.
///
/// Generates a single process vertex with the given number of tracks. Each
/// track's momentum and direction is drawn from uniform parameter
/// distributions. Position and time are always set to zero.
class ParametricParticleGenerator {
 public:
  struct Config {
    /// Number of particles.
    size_t numParticles = 1;
    /// Low, high for the transverse point of closest approach.
    std::array<double, 2> d0Range = {{0.0, 0.0}};
    /// Low, high for the z position at the point of closest approach.
    std::array<double, 2> z0Range = {{0.0, 0.0}};
    /// Low, high for the time at the point of closest approach.
    std::array<double, 2> t0Range = {{0.0, 0.0}};
    /// Low, high for the transverse angle.
    std::array<double, 2> phiRange = {{-M_PI, M_PI}};
    /// Low, high for pseudo-rapidity.
    std::array<double, 2> etaRange = {{-4.0, 4.0}};
    /// Low, high for transverse momentum.
    std::array<double, 2> ptRange = {
        {100 * Acts::UnitConstants::MeV, 10 * Acts::UnitConstants::GeV}};
    /// (Absolute) PDG particle number to identify the particle type.
    Acts::PdgParticle pdg = Acts::PdgParticle::eMuon;
    /// Randomize the charge and flip the PDG particle number sign accordingly.
    bool randomizeCharge = false;
  };

  ParametricParticleGenerator(const Config& cfg);

  /// Generate a single primary vertex with the given number of particles.
  SimParticleContainer operator()(RandomEngine& rng) const;

 private:
  Config m_cfg;
  // will be automatically set from PDG data tables
  double m_charge;
  double m_mass;
};

}  // namespace ActsExamples
