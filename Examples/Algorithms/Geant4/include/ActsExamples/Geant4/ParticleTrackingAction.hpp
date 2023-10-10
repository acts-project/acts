// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"

#include <memory>
#include <optional>
#include <string>

#include <G4Track.hh>
#include <G4UserTrackingAction.hh>

class G4Track;

namespace ActsExamples {

/// The G4UserTrackingAction that is called for every track in
/// the simulation process.
///
/// It records the initial and final particle state
class ParticleTrackingAction : public G4UserTrackingAction {
 public:
  struct Config {
    std::shared_ptr<EventStore> eventStore;

    bool keepParticlesWithoutHits = true;
  };

  /// Construct the stepping action
  ///
  /// @param cfg the configuration struct
  /// @param logger the ACTS logging instance
  ParticleTrackingAction(const Config& cfg,
                         std::unique_ptr<const Acts::Logger> logger =
                             Acts::getDefaultLogger("ParticleTrackingAction",
                                                    Acts::Logging::INFO));
  ~ParticleTrackingAction() override = default;

  /// Action before the track is processed in the
  /// the simulation, this will record the initial particle
  ///
  /// @param aTrack the current Geant4 track
  void PreUserTrackingAction(const G4Track* aTrack) final;

  /// Action after the track is processed in the
  /// the simulation, this will record the final particle
  ///
  /// @param aTrack the current Geant4 track
  void PostUserTrackingAction(const G4Track* aTrack) final;

 protected:
  Config m_cfg;

 private:
  /// Convert a G4Track to a SimParticle
  ///
  /// @param aTrack the current Geant4 track
  /// @param particleId the particle ID the particle will have
  SimParticle convert(const G4Track& aTrack, SimBarcode particleId) const;

  /// Make the particle id
  std::optional<SimBarcode> makeParticleId(G4int trackId, G4int parentId) const;

  /// Private access method to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// Private access method to the event store
  EventStore& eventStore() const { return *m_cfg.eventStore; }

  /// The looging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples
