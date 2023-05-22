// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Geant4/EventStoreRegistry.hpp"

#include <G4ParticleDefinition.hh>
#include <G4RunManager.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>
#include <globals.hh>

ActsExamples::ParticleTrackingAction::ParticleTrackingAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserTrackingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void ActsExamples::ParticleTrackingAction::PreUserTrackingAction(
    const G4Track* aTrack) {
  auto& eventData = EventStoreRegistry::eventData();

  auto particleId = makeParticleId(aTrack->GetTrackID(), aTrack->GetParentID());

  // There is already a warning printed in the makeParticleId function
  if (not particleId) {
    return;
  }

  auto [it, success] =
      eventData.particlesInitial.insert(convert(*aTrack, *particleId));

  // Only register particle at the initial state AND if there is no particle ID
  // collision
  if (success) {
    eventData.trackIdMapping[aTrack->GetTrackID()] = *particleId;
  } else {
    eventData.particleIdCollisionsInitial++;
    ACTS_WARNING("Particle ID collision with "
                 << *particleId
                 << " detected for initial particles. Skip particle");
  }
}

void ActsExamples::ParticleTrackingAction::PostUserTrackingAction(
    const G4Track* aTrack) {
  auto& eventData = EventStoreRegistry::eventData();

  // The initial particle maybe was not registered because a particle ID
  // collision
  if (eventData.trackIdMapping.find(aTrack->GetTrackID()) ==
      eventData.trackIdMapping.end()) {
    return;
  }

  const auto barcode = eventData.trackIdMapping.at(aTrack->GetTrackID());

  auto hasHits = eventData.particleHitCount.find(barcode) !=
                     eventData.particleHitCount.end() and
                 eventData.particleHitCount.at(barcode) > 0;

  if (not m_cfg.keepParticlesWithoutHits and not hasHits) {
    auto n = eventData.particlesInitial.erase(
        ActsExamples::SimParticle{barcode, Acts::PdgParticle::eInvalid});
    assert(n == 1);
    return;
  }

  auto particle = convert(*aTrack, barcode);
  auto [it, success] = eventData.particlesFinal.insert(particle);

  if (not success) {
    eventData.particleIdCollisionsFinal++;
    ACTS_WARNING("Particle ID collision with "
                 << particle.particleId()
                 << " detected for final particles. Skip particle");
  }
}

ActsExamples::SimParticle ActsExamples::ParticleTrackingAction::convert(
    const G4Track& aTrack, SimBarcode particleId) const {
  // Unit conversions G4->::ACTS
  constexpr double convertTime = Acts::UnitConstants::s / CLHEP::s;
  constexpr double convertLength = Acts::UnitConstants::mm / CLHEP::mm;
  constexpr double convertEnergy = Acts::UnitConstants::GeV / CLHEP::GeV;

  // Get all the information from the Track
  const G4ParticleDefinition* particleDef = aTrack.GetParticleDefinition();
  G4int pdg = particleDef->GetPDGEncoding();
  G4double charge = particleDef->GetPDGCharge();
  G4double mass = particleDef->GetPDGMass();
  G4ThreeVector pPosition = convertLength * aTrack.GetPosition();
  G4double pTime = convertTime * aTrack.GetGlobalTime();
  G4ThreeVector pDirection = aTrack.GetMomentumDirection();
  G4double p = convertEnergy * aTrack.GetKineticEnergy();

  // Now create the Particle
  ActsExamples::SimParticle aParticle(particleId, Acts::PdgParticle(pdg),
                                      charge, mass);
  aParticle.setPosition4(pPosition[0], pPosition[1], pPosition[2], pTime);
  aParticle.setDirection(pDirection[0], pDirection[1], pDirection[2]);
  aParticle.setAbsoluteMomentum(p);
  return aParticle;
}

std::optional<ActsExamples::SimBarcode>
ActsExamples::ParticleTrackingAction::makeParticleId(G4int trackId,
                                                     G4int parentId) const {
  auto& ed = EventStoreRegistry::eventData();

  // We already have this particle registered (it is one of the input particles
  // or we are making a final particle state)
  if (ed.trackIdMapping.find(trackId) != ed.trackIdMapping.end()) {
    return ed.trackIdMapping.at(trackId);
  }

  if (ed.trackIdMapping.find(parentId) == ed.trackIdMapping.end()) {
    ACTS_DEBUG("Parent particle " << parentId
                                  << " not registered, cannot build barcode");
    ed.parentIdNotFound++;
    return std::nullopt;
  }

  auto pid = ed.trackIdMapping.at(parentId).makeDescendant();

  auto key = EventStoreRegistry::State::BarcodeWithoutSubparticle::Zeros();
  key.set(0, pid.vertexPrimary())
      .set(1, pid.vertexSecondary())
      .set(2, pid.particle())
      .set(3, pid.generation());
  pid.setSubParticle(++ed.subparticleMap[key]);

  return pid;
}
