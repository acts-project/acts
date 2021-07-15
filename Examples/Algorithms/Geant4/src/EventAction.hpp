// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/MaterialInteractor.hpp"
#include "ActsExamples/EventData/SimHit.hpp"

#include <memory>

#include <G4UserEventAction.hh>
#include <globals.hh>

namespace ActsExamples::Geant4 {

class SteppingAction;

/// @class EventAction
///
/// @brief Writes out material track records
///
/// The EventAction class is the realization of the Geant4 class
/// G4UserEventAction and is writing out the collected RecordedMaterialTrack
/// entities needed for material mapping once per event.
///
class EventAction final : public G4UserEventAction {
 public:
  /// Static access method
  static EventAction* instance();

  /// Construct the action and ensure singleton usage.
  EventAction();
  ~EventAction() final override;

  /// Interface method for begin of the event
  /// @param event is the G4Event to be processed
  /// @note resets the material step action
  void BeginOfEventAction(const G4Event* event) final override;

  /// Interface method for end of event
  /// @param event is the G4Event to be processed
  /// @note this method is writing out the material track records
  void EndOfEventAction(const G4Event* event) final override;

  /// Clear the recorded data.
  void clear();

  /// Access the recorded material tracks.
  ///
  /// This only contains valid data after the end-of-event action has been
  /// executed.
  const std::vector<Acts::RecordedMaterialTrack>& materialTracks() const;

 private:
  /// Instance of the EventAction
  static EventAction* s_instance;

  /// The materialTrackWriter
  std::vector<Acts::RecordedMaterialTrack> m_materialTracks;
};

}  // namespace ActsExamples::Geant4
