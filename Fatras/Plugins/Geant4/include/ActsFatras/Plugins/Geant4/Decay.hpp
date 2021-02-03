// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Plugins/Geant4/PDGtoG4Converter.hpp"

#include <cmath>
#include <random>
#include <vector>

#include "G4ParticleDefinition.hh"

class G4RunManager;

namespace ActsFatras {

/// @brief This class handles the decay of a particle using the Geant4 classes
/// and provides the decay products. Additionally this class assigns a lifetime
/// to a particle. This lifetime is used to trigger the decay.
class Decay {
 public:
  using Scalar = Particle::Scalar;

  /// Deviation factor from generated proper time limit
  Scalar properTimeTolerance = 1e-2;
  
  /// Constructor
  Decay();

  /// @brief This method assigns a randomised lifetime to a particle based on
  /// its type.
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in] particle The particle that gets a lifetime assigned
  ///
  /// @return Proper time limit of the particle
  template <typename generator_t>
  Scalar generateProperTimeLimit(generator_t& generator,
                                 const Particle& particle) const;

  /// @brief This function tests if a decay should occur, triggers it whenever
  /// necessary and evaluates the decay products.
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] particle The particle that may decay
  ///
  /// @return Vector containing decay products
  template <typename generator_t>
  std::vector<Particle> run(generator_t&, Particle& particle) const;

  /// @brief This function evaluates the decay products of a given particle
  ///
  /// @param [in] parent The decaying particle
  ///
  /// @return Vector containing the decay products
  std::vector<Particle> decayParticle(const Particle& parent) const;

 private:
  /// @brief Convenience method assuring the existance of a G4RunManager
  ///
  /// @return Pointer to the run manager
  G4RunManager* initG4RunManager() const;

  mutable G4RunManager* m_g4RunManager;  ///< for dummy G4 initialization

  PDGtoG4Converter
      m_pdgToG4Conv;  ///< Handle for converting a PDG ID into a Geant4 particle
};

template <typename generator_t>
Particle::Scalar Decay::generateProperTimeLimit(
    generator_t& generator, const Particle& particle) const {
  // Get the particle properties
  const int pdgCode = particle.pdg();
  // Keep muons stable
  if (std::abs(pdgCode) == 13)
    return std::numeric_limits<Scalar>::infinity();

  // Get the Geant4 particle
  G4ParticleDefinition* pDef = m_pdgToG4Conv.getParticleDefinition(pdgCode);

  // Fast exit if the particle is stable
  if (!pDef || pDef->GetPDGStable()) {
    return std::numeric_limits<Scalar>::infinity();
  }

  // Get average lifetime
  constexpr Scalar convertTime = Acts::UnitConstants::mm / CLHEP::s;
  const Scalar tau = pDef->GetPDGLifeTime() * convertTime;
  // Sample & return the lifetime
  std::uniform_real_distribution<Scalar> uniformDistribution{0., 1.};

  return -tau * log(uniformDistribution(generator));
}

template <typename generator_t>
std::vector<Particle> Decay::run(generator_t&, Particle& particle) const {
  // Fast exit if particle is not alive
  if (!particle)
    return {};

  // Perform the decay
  std::vector<Particle> decayProducts = decayParticle(particle);

  // Kill the particle
  particle.setAbsoluteMomentum(Scalar(0));
  return decayProducts;
}
}  // namespace ActsFatras