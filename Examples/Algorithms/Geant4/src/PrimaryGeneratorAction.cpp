// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/PrimaryGeneratorAction.hpp"

#include <stdexcept>

#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4RandomDirection.hh>
#include <G4UnitsTable.hh>
#include <Randomize.hh>

namespace ActsExamples::Geant4 {

PrimaryGeneratorAction* PrimaryGeneratorAction::s_instance = nullptr;

PrimaryGeneratorAction* PrimaryGeneratorAction::instance() {
  return s_instance;
}

PrimaryGeneratorAction::PrimaryGeneratorAction(const Config& cfg)
    : G4VUserPrimaryGeneratorAction(),
      m_cfg(cfg),
      m_particleGun(std::make_unique<G4ParticleGun>(1)) {
  if (s_instance) {
    throw std::logic_error(
        "Attempted to duplicate the PrimaryGeneratorAction singleton");
  } else {
    s_instance = this;
  }

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle =
      particleTable->FindParticle(m_cfg.particleName);
  m_particleGun->SetParticleDefinition(particle);
  m_particleGun->SetParticleEnergy(m_cfg.energy);
  G4UnitDefinition::PrintUnitsTable();

  // set the random seeds
  CLHEP::HepRandom::getTheEngine()->setSeed(m_cfg.randomSeed1,
                                            m_cfg.randomSeed2);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  s_instance = nullptr;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  // this function is called at the begining of event
  G4double phi =
      m_cfg.phiRange.first +
      G4UniformRand() * (m_cfg.phiRange.second - m_cfg.phiRange.first);
  G4double theta = 0;
  if (m_cfg.samplingVariable == "theta") {
    G4double thetaMin = 2 * atan(exp(-m_cfg.etaRange.first));
    G4double thetaMax = 2 * atan(exp(-m_cfg.etaRange.second));
    theta = thetaMin + G4UniformRand() * (thetaMax - thetaMin);
  } else if (m_cfg.samplingVariable == "eta") {
    G4double eta =
        m_cfg.etaRange.first +
        G4UniformRand() * (m_cfg.etaRange.second - m_cfg.etaRange.first);
    theta = 2 * atan(exp(-eta));
  } else {
    throw std::invalid_argument("Unknow sampling variable");
  }
  // build a direction
  m_direction =
      G4ThreeVector(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
  m_position =
      G4ThreeVector(G4RandGauss::shoot(m_cfg.vertexPosX, m_cfg.vertexSigmaX),
                    G4RandGauss::shoot(m_cfg.vertexPosY, m_cfg.vertexSigmaY),
                    G4RandGauss::shoot(m_cfg.vertexPosZ, m_cfg.vertexSigmaZ));
  // set to the particle gun and
  m_particleGun->SetParticleMomentumDirection(m_direction);
  m_particleGun->SetParticlePosition(m_position);
  m_particleGun->GeneratePrimaryVertex(anEvent);
}

}  // namespace ActsExamples::Geant4