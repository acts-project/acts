// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/Plugins/Geant4/Geant4Decay.hpp"

#include "Acts/Definitions/Common.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"
#include "ActsFatras/Plugins/Geant4/G4DetectorConstruction.hpp"

#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4RunManager.hh"
#include "QGSP_BERT.hh"

ActsFatras::Geant4Decay::Geant4Decay() : m_g4RunManager(initG4RunManager()) {}

std::vector<ActsFatras::Particle> ActsFatras::Geant4Decay::decayParticle(
    const ActsFatras::Particle& parent) const {
  std::vector<Particle> children;

  // Find the particle type that will decay
  int pdgCode = parent.pdg();
  G4ParticleDefinition* pDef = m_pdgToG4Conv.getParticleDefinition(pdgCode);
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

G4RunManager* ActsFatras::Geant4Decay::initG4RunManager() const {
  // Test if there's already a G4RunManager
  if (G4RunManager::GetRunManager() == nullptr) {
    G4RunManager* runManager = new G4RunManager;

    // Initialise physics
    G4VUserPhysicsList* thePL = new QGSP_BERT;
    runManager->SetUserInitialization(thePL);

    // Build a dummy detector
    runManager->SetUserInitialization(new G4DetectorConstruction());

    // Initialise the G4RunManager itself
    runManager->Initialize();
    return runManager;
  } else {
    // Return the existing G4RunManager
    return G4RunManager::GetRunManager();
  }
}