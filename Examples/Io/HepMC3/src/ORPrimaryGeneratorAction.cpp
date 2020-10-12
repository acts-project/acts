// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ORPrimaryGeneratorAction.hpp"
#include <stdexcept>
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

ActsExamples::ORPrimaryGeneratorAction*
    ActsExamples::ORPrimaryGeneratorAction::s_instance
    = nullptr;

ActsExamples::ORPrimaryGeneratorAction::ORPrimaryGeneratorAction(G4int           randomSeed1,
    G4int           randomSeed2)
  : G4VUserPrimaryGeneratorAction(), m_particleGun(nullptr)
{
  // Configure the run
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    s_instance = this;
  }
  
  // Configure the gun
  G4int nofParticles = 1;
  m_particleGun       = std::make_unique<G4ParticleGun>(nofParticles);

  // Prepare the particle table
  m_particleTable = G4ParticleTable::GetParticleTable();
  
  // Set the random seeds
  CLHEP::HepRandom::getTheEngine()->setSeed(randomSeed1, randomSeed2);
}

ActsExamples::ORPrimaryGeneratorAction::~ORPrimaryGeneratorAction()
{
  s_instance = nullptr;
}

ActsExamples::ORPrimaryGeneratorAction*
ActsExamples::ORPrimaryGeneratorAction::instance()
{
  // Static acces function via G4RunManager
  return s_instance;
}

void ActsExamples::ORPrimaryGeneratorAction::prepareParticleGun(G4int pdg,
    G4double        momentum,
	 G4ThreeVector pos, G4ThreeVector dir)
{
	// Particle type
  G4ParticleDefinition* particle = m_particleTable->FindParticle(pdg);
  m_particleGun->SetParticleDefinition(particle);
  // Particle properties
  m_particleGun->SetParticlePosition(pos); 
  m_particleGun->SetParticleMomentum(momentum);
  m_particleGun->SetParticleMomentumDirection(dir);
}

void
ActsExamples::ORPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	// Produce the event
  m_particleGun->GeneratePrimaryVertex(anEvent);
}
