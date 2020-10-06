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
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "G4UserEventAction.hh"
#include "globals.hh"

namespace ActsExamples {
	
	struct Collection
	{
		//~ std::map<int, std::vector<ParticleRecord>> particles;
		int pdg;
		double momentum, phi, theta;
	};
	
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
	
	Collection
	//~ processTracks(const int pdg, const double momentum, const double phi, const double theta) const
	processTracks() const
	{
		// TODO
		return m_processTracks;
	}
	
	//~ Collection	
	//~ outcomingParticles(const int pdg, const double momentum, const double phi, const double theta) const
	//~ {
		//~ Collection c;
		//~ c.particles = m_particles;
		//~ c.pdg = pdg;
		//~ c.momentum = momentum;
		//~ c.phi = phi;
		//~ c.theta = theta;
		//~ return c;
	//~ }
    
	private:
	/// Instance of the EventAction
	static OREventAction* s_instance;

	Collection m_processTracks;	
};
}  // namespace ActsExamples