// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/ParticleTrackingAction.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/UnitConversion.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <cassert>
#include <ostream>
#include <unordered_map>
#include <utility>

#include <G4ParticleDefinition.hh>
#include <G4RunManager.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>

namespace ActsExamples::Geant4 {

ParticleTrackingAction::ParticleTrackingAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserTrackingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void ParticleTrackingAction::PreUserTrackingAction(const G4Track* trackPtr) {
  assert(trackPtr != nullptr);
  const G4Track& track = *trackPtr;

  // If this is not the case, there are unhandled cases of particle stopping in
  // the SensitiveSteppingAction
  // TODO We could also merge the remaining hits to a hit here, but it would be
  // nicer to investigate, if we can handle all particle stop conditions in the
  // SensitiveSteppingAction... This seems to happen O(1) times in a ttbar
  // event, so seems not to be too problematic
  if (!eventStore().hitBuffer.empty()) {
    eventStore().hitBuffer.clear();
    ACTS_WARNING("Hit buffer not empty after track");
  }

  const std::optional<SimBarcode> barcode =
      makeParticleId(track.GetTrackID(), track.GetParentID());

  // There is already a warning printed in the makeParticleId function if this
  // indicates a failure
  if (!barcode.has_value()) {
    return;
  }

  const SimParticleState fatrasParticle = convert(track, *barcode);
  const SimParticle particle(fatrasParticle, fatrasParticle);
  const auto [it, success] = eventStore().particlesInitial.insert(particle);

  // Only register particle at the initial state AND if there is no particle ID
  // collision
  if (success) {
    eventStore().trackIdMapping[track.GetTrackID()] = particle.particleId();
  } else {
    eventStore().particleIdCollisionsInitial++;
    ACTS_WARNING("Particle ID collision with "
                 << particle.particleId()
                 << " detected for initial particles. Skip particle");
  }
}

void ParticleTrackingAction::PostUserTrackingAction(const G4Track* trackPtr) {
  assert(trackPtr != nullptr);
  const G4Track& track = *trackPtr;

  // The initial particle maybe was not registered because of a particle ID
  // collision
  if (!eventStore().trackIdMapping.contains(track.GetTrackID())) {
    ACTS_WARNING("Particle ID for track ID " << track.GetTrackID()
                                             << " not registered. Skip");
    return;
  }

  const SimBarcode barcode = eventStore().trackIdMapping.at(track.GetTrackID());

  const bool hasHits = eventStore().particleHitCount.contains(barcode) &&
                       eventStore().particleHitCount.at(barcode) > 0;
  if (!m_cfg.keepParticlesWithoutHits && !hasHits) {
    [[maybe_unused]] const std::size_t n =
        eventStore().particlesSimulated.erase(
            SimParticle(barcode, Acts::PdgParticle::eInvalid));
    assert(n == 1);
    return;
  }

  const auto particleIt = eventStore().particlesInitial.find(barcode);
  if (particleIt == eventStore().particlesInitial.end()) {
    ACTS_WARNING("Particle ID " << barcode
                                << " not found in initial particles");
    return;
  }
  SimParticle particle = *particleIt;
  particle.finalState() = convert(track, barcode);

  const auto [it, success] = eventStore().particlesSimulated.insert(particle);

  if (!success) {
    eventStore().particleIdCollisionsFinal++;
    ACTS_WARNING("Particle ID collision with "
                 << particle.particleId()
                 << " detected for final particles. Skip particle");
  }
}

SimParticleState ParticleTrackingAction::convert(const G4Track& track,
                                                 SimBarcode particleId) const {
  // Get all the information from the Track
  const G4ParticleDefinition* particleDef = track.GetParticleDefinition();
  const G4int pdg = particleDef->GetPDGEncoding();
  const G4double charge = particleDef->GetPDGCharge();
  const G4double mass = convertEnergyToActs * particleDef->GetPDGMass();
  const G4ThreeVector pPosition = convertLengthToActs * track.GetPosition();
  const G4double pTime = convertTimeToActs * track.GetGlobalTime();
  const G4ThreeVector pDirection = track.GetMomentumDirection();
  const G4double p = convertEnergyToActs * track.GetKineticEnergy();

  std::uint32_t numberOfHits = 0;
  if (const auto it = eventStore().particleHitCount.find(particleId);
      it != eventStore().particleHitCount.end()) {
    numberOfHits = it->second;
  }

  ActsFatras::ParticleOutcome particleOutcome =
      ActsFatras::ParticleOutcome::Alive;
  if (const auto it = eventStore().particleOutcome.find(particleId);
      it != eventStore().particleOutcome.end()) {
    particleOutcome = it->second;
  }

  // Now create the Particle
  SimParticleState aParticle(particleId, Acts::PdgParticle{pdg}, charge, mass);
  aParticle.setPosition4(pPosition[0], pPosition[1], pPosition[2], pTime);
  aParticle.setDirection(pDirection[0], pDirection[1], pDirection[2]);
  aParticle.setAbsoluteMomentum(p);
  aParticle.setNumberOfHits(numberOfHits);
  aParticle.setOutcome(particleOutcome);
  return aParticle;
}

std::optional<SimBarcode> ParticleTrackingAction::makeParticleId(
    G4int trackId, G4int parentId) const {
  // We already have this particle registered (it is one of the input particles
  // or we are making a final particle state)
  if (eventStore().trackIdMapping.contains(trackId)) {
    return std::nullopt;
  }

  if (!eventStore().trackIdMapping.contains(parentId)) {
    ACTS_DEBUG("Parent particle " << parentId
                                  << " not registered, cannot build barcode");
    eventStore().parentIdNotFound++;
    return std::nullopt;
  }

  SimBarcode pid = eventStore().trackIdMapping.at(parentId).makeDescendant();
  const SimBarcode key = pid.withoutSubparticle();
  ++eventStore().subparticleMap[key];
  pid = pid.withSubParticle(eventStore().subparticleMap[key]);

  return pid;
}

}  // namespace ActsExamples::Geant4
