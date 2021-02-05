// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"
#include "ActsFatras/Geant4/DummyDetectorConstruction.hpp"
#include "ActsFatras/Geant4/PDGtoG4Converter.hpp"

#include <cmath>
#include <random>
#include <vector>

#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4ParticleDefinition.hh"

class G4RunManager;

namespace ActsFatras {

/// Handle particle decays using the Geant4 decay models.
class Geant4Decay {
 public:
  using Scalar = Particle::Scalar;

  /// Constructor
  Geant4Decay() : : m_g4RunManager(ensureGeant4RunManager()) {}

  /// Generate the particle life time.
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in] particle The particle that gets a lifetime assigned
  ///
  /// @return Proper time limit of the particle
  template <typename generator_t>
  Scalar generateProperTimeLimit(generator_t& generator,
                                 const Particle& particle) const;

  /// Decay the particle and create the decay products.
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] particle The particle that may decay
  ///
  /// @return Vector containing decay products
  template <typename generator_t>
  std::vector<Particle> run(generator_t&, Particle& particle) const;

 private:
  mutable G4RunManager* m_g4RunManager;  ///< for dummy G4 initialization

  PDGtoG4Converter
      m_pdgToG4Conv;  ///< Handle for converting a PDG ID into a Geant4 particle
};

template <typename generator_t>
Particle::Scalar Geant4Decay::generateProperTimeLimit(
    generator_t& generator, const Particle& particle) const {
  // Get the particle properties
  const Acts::PdgParticle pdgCode = particle.pdg();
  // Keep muons stable
  if (makeAbsolutePdgParticle(pdgCode) == Acts::PdgParticle::eMuon)
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

  return -tau * std::log(uniformDistribution(generator));
}

template <typename generator_t>
std::vector<Particle> Geant4Decay::run(generator_t&, Particle& particle) const {
  std::vector<Particle> children;

  // Find the particle type that will decay
  G4ParticleDefinition* pDef =
      m_pdgToG4Conv.getParticleDefinition(parent.pdg());
  if (!pDef) {
    return children;
  }
  // Get the particles decay table
  G4DecayTable* dt = pDef->GetDecayTable();
  if (!dt) {
    return children;
  }
  // Select a decay channel
  G4VDecayChannel* channel = dt->SelectADecayChannel();
  if (!channel) {
    return children;
  }
  // Get the decay products from the selected channel
  G4DecayProducts* products = channel->DecayIt();
  if (!products) {
    return children;
  }

  // Boost the decay products using the parents four-momentum
  const Particle::Vector4 mom4 = parent.fourMomentum();
  products->Boost(mom4[Acts::eMom0] / mom4[Acts::eEnergy],
                  mom4[Acts::eMom1] / mom4[Acts::eEnergy],
                  mom4[Acts::eMom2] / mom4[Acts::eEnergy]);

  G4int nProducts = products->entries();
  for (G4int i = 0; i < nProducts; i++) {
    G4DynamicParticle* prod = products->PopProducts();
    if (!prod) {
      continue;
    }

    // Convert the decay product from Geant4 to Acts
    const G4ThreeVector& mom = prod->GetMomentum();
    constexpr Scalar convertEnergy = Acts::UnitConstants::GeV / CLHEP::GeV;
    Acts::Vector3 amgMom(mom.x(), mom.y(), mom.z());
    amgMom *= convertEnergy;
    const int32_t pdg = prod->GetPDGcode();

    Particle childParticle(Barcode(), static_cast<Acts::PdgParticle>(pdg));
    childParticle.setPosition4(parent.fourPosition())
        .setAbsoluteMomentum(amgMom.norm())
        .setDirection(amgMom)
        .setProcess(ProcessType::eDecay);

    // Store the particle
    children.push_back(std::move(childParticle));
  }
  return children;
}
}  // namespace ActsFatras