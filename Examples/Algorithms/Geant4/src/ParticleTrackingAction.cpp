// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/MultiIndex.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <cassert>
#include <cstddef>
#include <ostream>
#include <unordered_map>
#include <utility>

#include <G4ParticleDefinition.hh>
#include <G4RunManager.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>

ActsExamples::ParticleTrackingAction::ParticleTrackingAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserTrackingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void ActsExamples::ParticleTrackingAction::PreUserTrackingAction(
    const G4Track* aTrack) {
  // If this is not the case, there are unhandled cases of particle stopping in
  // the SensitiveSteppingAction
  // TODO We could also merge the remaining hits to a hit here, but it would be
  // nicer to investigate, if we can handle all particle stop conditions in the
  // SensitiveSteppingAction... This seems to happen O(1) times in a ttbar
  // event, so seems not to be too problematic
  if (not eventStore().hitBuffer.empty()) {
    eventStore().hitBuffer.clear();
    ACTS_WARNING("Hit buffer not empty after track");
  }

  auto particleId = makeParticleId(aTrack->GetTrackID(), aTrack->GetParentID());

  // There is already a warning printed in the makeParticleId function
  if (not particleId) {
    return;
  }

  auto [it, success] =
      eventStore().particlesInitial.insert(convert(*aTrack, *particleId));

  // Only register particle at the initial state AND if there is no particle ID
  // collision
  if (success) {
    eventStore().trackIdMapping[aTrack->GetTrackID()] = *particleId;
  } else {
    eventStore().particleIdCollisionsInitial++;
    ACTS_WARNING("Particle ID collision with "
                 << *particleId
                 << " detected for initial particles. Skip particle");
  }
}

void ActsExamples::ParticleTrackingAction::PostUserTrackingAction(
    const G4Track* aTrack) {
  // The initial particle maybe was not registered because a particle ID
  // collision
  if (eventStore().trackIdMapping.find(aTrack->GetTrackID()) ==
      eventStore().trackIdMapping.end()) {
    return;
  }

  const auto barcode = eventStore().trackIdMapping.at(aTrack->GetTrackID());

  auto hasHits = eventStore().particleHitCount.find(barcode) !=
                     eventStore().particleHitCount.end() and
                 eventStore().particleHitCount.at(barcode) > 0;

  if (not m_cfg.keepParticlesWithoutHits and not hasHits) {
    [[maybe_unused]] auto n = eventStore().particlesInitial.erase(
        ActsExamples::SimParticle{barcode, Acts::PdgParticle::eInvalid});
    assert(n == 1);
    return;
  }

  auto particle = convert(*aTrack, barcode);
  auto [it, success] = eventStore().particlesFinal.insert(particle);

  if (not success) {
    eventStore().particleIdCollisionsFinal++;
    ACTS_WARNING("Particle ID collision with "
                 << particle.particleId()
                 << " detected for final particles. Skip particle");
  }
}

ActsExamples::SimParticle ActsExamples::ParticleTrackingAction::convert(
    const G4Track& aTrack, SimBarcode particleId) const {
  // Unit conversions G4->::ACTS
  constexpr double convertTime = Acts::UnitConstants::ns / CLHEP::ns;
  constexpr double convertLength = Acts::UnitConstants::mm / CLHEP::mm;
  constexpr double convertEnergy = Acts::UnitConstants::GeV / CLHEP::GeV;

  // Get all the information from the Track
  const G4ParticleDefinition* particleDef = aTrack.GetParticleDefinition();
  G4int pdg = particleDef->GetPDGEncoding();
  G4double charge = particleDef->GetPDGCharge();
  G4double mass = convertEnergy * particleDef->GetPDGMass();
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
  // We already have this particle registered (it is one of the input particles
  // or we are making a final particle state)
  if (eventStore().trackIdMapping.find(trackId) !=
      eventStore().trackIdMapping.end()) {
    return eventStore().trackIdMapping.at(trackId);
  }

  if (eventStore().trackIdMapping.find(parentId) ==
      eventStore().trackIdMapping.end()) {
    ACTS_DEBUG("Parent particle " << parentId
                                  << " not registered, cannot build barcode");
    eventStore().parentIdNotFound++;
    return std::nullopt;
  }

  auto pid = eventStore().trackIdMapping.at(parentId).makeDescendant();

  auto key = EventStore::BarcodeWithoutSubparticle::Zeros();
  key.set(0, pid.vertexPrimary())
      .set(1, pid.vertexSecondary())
      .set(2, pid.particle())
      .set(3, pid.generation());
  pid.setSubParticle(++eventStore().subparticleMap[key]);

  return pid;
}
