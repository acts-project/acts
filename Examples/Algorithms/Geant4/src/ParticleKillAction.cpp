// This file is part of the Acts project.
//
// Copyright (C) 2021-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/ParticleKillAction.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Geant4/Geant4ParticleStatus.hpp"

#include <ostream>
#include <utility>

#include <G4RunManager.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4Track.hh>
#include <G4UnitsTable.hh>
#include <G4VPhysicalVolume.hh>

ActsExamples::ParticleKillAction::ParticleKillAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserSteppingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void ActsExamples::ParticleKillAction::UserSteppingAction(const G4Step* step) {
  constexpr double convertLength = Acts::UnitConstants::mm / CLHEP::mm;
  constexpr double convertTime = Acts::UnitConstants::ns / CLHEP::ns;

  G4Track* track = step->GetTrack();

  const auto pos = convertLength * track->GetPosition();
  const auto time = convertTime * track->GetGlobalTime();
  const bool isSecondary =
      track->GetDynamicParticle()->GetPrimaryParticle() == nullptr;

  const bool outOfVolume = m_cfg.volume && !m_cfg.volume->inside(Acts::Vector3{
                                               pos.x(), pos.y(), pos.z()});
  const bool outOfTime = time > m_cfg.maxTime;
  const bool invalidSecondary = m_cfg.secondaries && isSecondary;

  if (outOfVolume || outOfTime || invalidSecondary) {
    ACTS_DEBUG("Kill track with internal track ID "
               << track->GetTrackID() << " at " << pos << " and global time "
               << time / Acts::UnitConstants::ns << "ns and isSecondary "
               << isSecondary);
    track->SetTrackStatus(G4TrackStatus::fStopAndKill);
  }

  // store the status of the particle
  // check if we have a particle assigned to track
  if (eventStore().trackIdMapping.find(track->GetTrackID()) ==
      eventStore().trackIdMapping.end()) {
    return;
  }

  // set the outcome of the particle
  const auto particleId = eventStore().trackIdMapping.at(track->GetTrackID());
  if (outOfVolume) {
    eventStore().particleStatus[particleId] =
        static_cast<std::uint32_t>(Geant4ParticleStatus::KilledVolume);
  } else if (outOfTime) {
    eventStore().particleStatus[particleId] =
        static_cast<std::uint32_t>(Geant4ParticleStatus::KilledTime);
  } else if (invalidSecondary) {
    eventStore().particleStatus[particleId] =
        static_cast<std::uint32_t>(Geant4ParticleStatus::KilledSecondary);
  } else if (track->GetTrackStatus() == fStopAndKill) {
    eventStore().particleStatus[particleId] =
        static_cast<std::uint32_t>(Geant4ParticleStatus::KilledAbsorbed);
  }
}
