// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/SensitiveSteppingAction.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Geant4/EventStoreRegistry.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"

#include <G4RunManager.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>
#include <G4VPhysicalVolume.hh>

ActsExamples::SensitiveSteppingAction::SensitiveSteppingAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserSteppingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void ActsExamples::SensitiveSteppingAction::UserSteppingAction(
    const G4Step* step) {
  // Unit conversions G4->::ACTS
  constexpr double convertTime = Acts::UnitConstants::s / CLHEP::s;
  constexpr double convertLength = Acts::UnitConstants::mm / CLHEP::mm;
  constexpr double convertEnergy = Acts::UnitConstants::GeV / CLHEP::GeV;

  // Get G4 pointers
  const G4Track* track = step->GetTrack();
  const G4PrimaryParticle* primaryParticle =
      track->GetDynamicParticle()->GetPrimaryParticle();
  const G4VPhysicalVolume* volume = track->GetVolume();

  // Retrieve the event data registry
  auto& eventData = EventStoreRegistry::eventData();

  // This is not the case if we have a particle-ID collision
  if (eventData.trackIdMapping.find(track->GetTrackID()) ==
      eventData.trackIdMapping.end()) {
    ACTS_WARNING("Probably we found a particle ID collision");
    return;
  }

  // Get the particle ID
  const auto particleID = eventData.trackIdMapping.at(track->GetTrackID());

  // Merge last hit if necessary. This is the case, if we are now in a different
  // volume then the volume in which the hit merger has been created, or if we
  // are tracking a different particle now
  auto& hitMerger = eventData.hitMerger;
  const bool resetHitMerger = (hitMerger->geant4Volume() != volume) or
                              (hitMerger->particleID() != particleID);

  if (hitMerger && resetHitMerger) {
    auto& hitCount = eventData.particleHitCount;

    // Get the hit count of the particle
    const auto index = hitCount.find(particleID) != hitCount.end()
                           ? hitCount.at(particleID)
                           : 0;

    // Merge the hit with a given hit index
    eventData.hitMergerSumHits += hitMerger->hits().size();
    const auto hit = hitMerger->mergeHit(index);

    // Kill hit merger after we extracted the merged hit
    hitMerger.reset();

    if (hit.depositedEnergy() > m_cfg.hitMergerEnergyThreshold) {
      // Increase counter (starts at 1 because of ++)
      ++hitCount[particleID];

      eventData.hits.push_back(hit);
    }
  }

  // Bail out if charged & configured to do so
  G4double absCharge = std::abs(track->GetParticleDefinition()->GetPDGCharge());
  if (not m_cfg.charged and absCharge > 0.) {
    return;
  }

  // Bail out if neutral & configured to do so
  if (not m_cfg.neutral and absCharge == 0.) {
    return;
  }

  // Bail out if it is a primary & configured to be ignored
  if (not m_cfg.primary and primaryParticle != nullptr) {
    return;
  }

  // Bail out if it is a secondary & configured to be ignored
  if (not m_cfg.secondary and primaryParticle == nullptr) {
    return;
  }

  // check if the volume has the sensitive string name
  const std::string_view volumeName(volume->GetName());
  if (volumeName.find(SensitiveSurfaceMapper::mappingPrefix) ==
      std::string_view::npos) {
    return;
  }

  // Cast out the GeometryIdentifier
  std::string volumeNameStr(volumeName);
  volumeNameStr.erase(0, SensitiveSurfaceMapper::mappingPrefix.size());

  const Acts::GeometryIdentifier geoID(std::stoul(volumeNameStr));

  ACTS_VERBOSE("Step of " << particleID << " in senstive volume "
                          << volumeName);

  // Get PreStepPoint and PostStepPoint
  const G4StepPoint* preStepPoint = step->GetPreStepPoint();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();

  G4ThreeVector preStepPosition = convertLength * preStepPoint->GetPosition();
  G4double preStepTime = convertTime * preStepPoint->GetGlobalTime();
  G4ThreeVector postStepPosition = convertLength * postStepPoint->GetPosition();
  G4double postStepTime = convertTime * postStepPoint->GetGlobalTime();

  G4ThreeVector preStepMomentum = convertEnergy * preStepPoint->GetMomentum();
  G4double preStepEnergy = convertEnergy * preStepPoint->GetTotalEnergy();
  G4ThreeVector postStepMomentum = convertEnergy * postStepPoint->GetMomentum();
  G4double postStepEnergy = convertEnergy * postStepPoint->GetTotalEnergy();

  Acts::ActsScalar hX = 0.5 * (preStepPosition[0] + postStepPosition[0]);
  Acts::ActsScalar hY = 0.5 * (preStepPosition[1] + postStepPosition[1]);
  Acts::ActsScalar hZ = 0.5 * (preStepPosition[2] + postStepPosition[2]);
  Acts::ActsScalar hT = 0.5 * (preStepTime + postStepTime);

  Acts::ActsScalar mXpre = preStepMomentum[0];
  Acts::ActsScalar mYpre = preStepMomentum[1];
  Acts::ActsScalar mZpre = preStepMomentum[2];
  Acts::ActsScalar mEpre = preStepEnergy;
  Acts::ActsScalar mXpost = postStepMomentum[0];
  Acts::ActsScalar mYpost = postStepMomentum[1];
  Acts::ActsScalar mZpost = postStepMomentum[2];
  Acts::ActsScalar mEpost = postStepEnergy;

  Acts::Vector4 particlePosition(hX, hY, hZ, hT);
  Acts::Vector4 beforeMomentum(mXpre, mYpre, mZpre, mEpre);
  Acts::Vector4 afterMomentum(mXpost, mYpost, mZpost, mEpost);

  // If we have no hit merger, create one
  if (not eventData.hitMerger) {
    eventData.hitMerger = HitMerger(particleID, geoID, volume);
  }

  eventData.hitMerger->addHit(ActsFatras::Hit(
      geoID, particleID, particlePosition, beforeMomentum, afterMomentum));
}
