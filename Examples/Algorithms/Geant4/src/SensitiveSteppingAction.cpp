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

namespace {
ActsFatras::Hit hitFromStep(const G4StepPoint* preStepPoint,
                            const G4StepPoint* postStepPoint,
                            ActsFatras::Barcode particleId,
                            Acts::GeometryIdentifier geoId, int32_t index) {
  static constexpr double convertTime = Acts::UnitConstants::s / CLHEP::s;
  static constexpr double convertLength = Acts::UnitConstants::mm / CLHEP::mm;
  static constexpr double convertEnergy = Acts::UnitConstants::GeV / CLHEP::GeV;

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

  return ActsFatras::Hit(geoId, particleId, particlePosition, beforeMomentum,
                         afterMomentum, index);
}
}  // namespace

ActsExamples::SensitiveSteppingAction::SensitiveSteppingAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserSteppingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void ActsExamples::SensitiveSteppingAction::UserSteppingAction(
    const G4Step* step) {
  // Unit conversions G4->::ACTS

  // Retrieve the event data registry
  auto& eventData = EventStoreRegistry::eventData();

  // The particle after the step
  G4Track* track = step->GetTrack();
  G4PrimaryParticle* primaryParticle =
      track->GetDynamicParticle()->GetPrimaryParticle();

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

  // Get the physical volume & check if it has the sensitive string name
  std::string_view volumeName = track->GetVolume()->GetName();

  if (volumeName.find(SensitiveSurfaceMapper::mappingPrefix) ==
      std::string_view::npos) {
    return;
  }

  ACTS_VERBOSE("Step in senstive volume " << volumeName);

  // Cast out the GeometryIdentifier
  volumeName = volumeName.substr(SensitiveSurfaceMapper::mappingPrefix.size(),
                                 std::string_view::npos);
  char* end = nullptr;
  const Acts::GeometryIdentifier geoId(
      std::strtoul(volumeName.data(), &end, 10));

  // This is not the case if we have a particle-ID collision
  if (eventData.trackIdMapping.find(track->GetTrackID()) ==
      eventData.trackIdMapping.end()) {
    return;
  }

  const auto particleId = eventData.trackIdMapping.at(track->GetTrackID());

  // Get PreStepPoint and PostStepPoint
  const G4StepPoint* preStepPoint = step->GetPreStepPoint();
  const G4StepPoint* postStepPoint = step->GetPostStepPoint();

  // Set particle hit count to zero, so we have this entry in the map later
  if (eventData.particleHitCount.find(particleId) ==
      eventData.particleHitCount.end()) {
    eventData.particleHitCount[particleId] = 0;
  }

  // Case A: The step starts at the entry of the volume and ends at the exit.
  // Add hit to collection.
  if (preStepPoint->GetStepStatus() == fGeomBoundary and
      postStepPoint->GetStepStatus() == fGeomBoundary) {
    const auto hit = hitFromStep(preStepPoint, postStepPoint, particleId, geoId,
                                 eventData.particleHitCount.at(particleId));

    if (hit.depositedEnergy() > m_cfg.hitMergerEnergyThreshold) {
      ++eventData.particleHitCount[particleId];
      eventData.hits.push_back(hit);
      eventData.numberGeantSteps += 1;
      eventData.maxStepsForHit = std::max(eventData.maxStepsForHit, 1ul);
    }

    return;
  }

  // Case B: The step doesn't end at the exit of the volume. Add the hit to the
  // hit buffer.
  if (postStepPoint->GetStepStatus() != fGeomBoundary) {
    eventData.hitBuffer.push_back(
        hitFromStep(preStepPoint, postStepPoint, particleId, geoId, -1));
    return;
  }

  // Case C: The step ends at the exit of the volume. Add hit to hit buffer, and
  // then combine hit buffer
  if (postStepPoint->GetStepStatus() == fGeomBoundary) {
    auto& buffer = eventData.hitBuffer;
    buffer.push_back(
        hitFromStep(preStepPoint, postStepPoint, particleId, geoId, -1));

    const auto pos4 =
        0.5 * (buffer.front().fourPosition() + buffer.back().fourPosition());

    const ActsFatras::Hit hit(geoId, particleId, pos4,
                              buffer.front().momentum4Before(),
                              buffer.back().momentum4After(),
                              eventData.particleHitCount.at(particleId));

    if (hit.depositedEnergy() > m_cfg.hitMergerEnergyThreshold) {
      ++eventData.particleHitCount[particleId];
      eventData.hits.push_back(hit);
      eventData.numberGeantSteps += buffer.size();
      eventData.maxStepsForHit =
          std::max(eventData.maxStepsForHit, buffer.size());
    }

    buffer.clear();
    return;
  }

  assert(false && "should never reach this");
}
