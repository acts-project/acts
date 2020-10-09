// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <map>
#include <string>
#include "G4UserEventAction.hh"
#include "globals.hh"

#include <HepMC3/GenEvent.h>

namespace ActsExamples {
	
/// @class EventAction
///
/// @brief Writes out material track records
///
/// The EventAction class is the realization of the Geant4 class
/// G4UserEventAction and is writing out the collected RecordedMaterialTrack
/// entities needed for material mapping once per event.
///
class OREventAction final : public G4UserEventAction {
 public:
  /// Static access method
  static OREventAction* instance();

  /// Construct the action and ensure singleton usage.
  OREventAction();
  ~OREventAction() final override;

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
  
  std::shared_ptr<HepMC3::GenEvent>
  event()
  {
	  return m_event;
  }
    
	private:
	/// Instance of the EventAction
	static OREventAction* s_instance;
    
    std::shared_ptr<HepMC3::GenEvent> m_event = nullptr;
};
}  // namespace ActsExamples