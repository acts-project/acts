// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/Physics/Decay/PDGtoG4Converter.hpp"
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
  /// Constructor
  Decay();

  /// @brief This method assigns a randomised lifetime to a particle based on
  /// its type.
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in, out] isp The particle that gets a lifetime assigned
  template <typename generator_t>
  void lifeTime(generator_t& generator, Particle& isp) const;

  /// @brief This function tests if a decay should occur, triggers it whenever
  /// necessary and evaluates the decay products.
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in] slab The material slab
  /// @param [in, out] isp The particle that may decay
  ///
  /// @return Vector containing decay products
  template <typename generator_t>
  std::vector<Particle> operator()(generator_t& generator,
                                   const Acts::MaterialSlab& /*slab*/,
                                   Particle& isp) const;

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

}  // namespace ActsFatras

template <typename generator_t>
void ActsFatras::Decay::lifeTime(generator_t& generator,
                                 ActsFatras::Particle& isp) const {
  // Get the particle properties
  const int pdgCode = isp.pdg();
  // Keep muons stable
  if (std::abs(pdgCode) == 13)
    return;

  // Get the Geant4 particle
  G4ParticleDefinition* pDef = m_pdgToG4Conv.getParticleDefinition(pdgCode);

  // Fast exit if the particle is stable
  if (!pDef || pDef->GetPDGStable()) {
    return;
  }

  // Get average lifetime
  constexpr double convertTime = Acts::UnitConstants::mm / CLHEP::s;
  const double tau = pDef->GetPDGLifeTime() * convertTime;
  // Sample the lifetime
  std::uniform_real_distribution<double> uniformDistribution{0., 1.};
  const double lifeTime = -tau * log(uniformDistribution(generator));

  // Assign the lifetime
  isp.setLifetimeLimit(lifeTime);
}

template <typename generator_t>
std::vector<ActsFatras::Particle> ActsFatras::Decay::operator()(
    generator_t& generator, const Acts::MaterialSlab& /*slab*/,
    ActsFatras::Particle& isp) const {
  // Test if decay condition is fulfilled
  if (isp.pathLimitTime() < isp.time())
    return {};

  // perform the decay
  std::vector<Particle> decayProducts = decayParticle(isp);

  // Assign a lifetime to decay products
  for (Particle& decayProduct : decayProducts) {
    lifeTime(generator, decayProduct);
  }

  // Kill the particle
  isp.setAbsMomentum(Particle::Scalar(0));

  return decayProducts;
}