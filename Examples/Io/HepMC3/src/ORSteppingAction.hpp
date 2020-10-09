// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include "G4UserSteppingAction.hh"
#include "globals.hh"

#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

namespace ActsExamples {

	/// @brief Data container of a particle
	// TODO: maybe replace this by an Acts particle
	struct ParticleRecord
	{
		std::array<double, 3> position, momentum, vertex;
		int pdg;
		double energy, mass, globalTime;
		int charge, trackid, parentid;
		std::string volume, process;
	};
	
  /// @class ORSteppingAction
  ///
  /// @brief Collects the particle
  class ORSteppingAction : public G4UserSteppingAction
  {
  public:
    ORSteppingAction();
    ~ORSteppingAction() override;

    /// Static access method to the instance
    static ORSteppingAction*
    instance();

    /// @brief Interface Method doing the step
    /// @note it creates and collects the MaterialInteraction entities
    /// @param step is the Geant4 step of the particle
    void
    UserSteppingAction(const G4Step* step) final override;

    /// Interface reset method
    /// @note it clears the collected step vector
    void
    clear();
    
    const std::map<int, std::vector<ParticleRecord>> processSteps() const {return m_particles;}

unsigned int& counter() {return testCounter;} 
unsigned int testCounter = 0;
  private:	
    /// Instance of the SteppingAction
    static ORSteppingAction* s_instance;
        
    std::shared_ptr<HepMC3::GenVertex> m_previousVertex = nullptr;    
        
	std::map<int, std::vector<ParticleRecord>> m_particles;	
  };
}  // namespace ActsExamples