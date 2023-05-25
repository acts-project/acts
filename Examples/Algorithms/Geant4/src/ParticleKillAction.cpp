// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/ParticleKillAction.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Geant4/EventStoreRegistry.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"

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

  const bool outOfVolume =
      m_cfg.volume and
      not m_cfg.volume->inside(Acts::Vector3{pos.x(), pos.y(), pos.z()});
  const bool outOfTime = time > m_cfg.maxTime;

  if (outOfVolume or outOfTime) {
    ACTS_DEBUG("Kill track with internal track ID "
               << track->GetTrackID() << " at " << pos << " and global time "
               << time / Acts::UnitConstants::ns << "ns");
    track->SetTrackStatus(G4TrackStatus::fStopAndKill);
  }
}
