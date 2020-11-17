// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <unordered_map>

class G4ParticleDefinition;

namespace ActsFatras
{
	/// @brief This class converts a PDG ID into a corresponding Geant4 particle.
  class PDGtoG4Converter
  {
  public:
    /// Constructor
    PDGtoG4Converter();

/// @brief Converts a PDG ID into the corresponding Geant4 particle
///
/// @param [in] pdgCode The PDG ID
///
/// @return Pointer to the Geant4 particle
G4ParticleDefinition* getParticleDefinition( int pdgCode) const;

  private:
    /// @brief This method fills the internal storage with Geant4 particles and their PDG IDs for later lookup
    void fillPredefinedParticles();
    
/// @brief This method adds a certain Particle to the internal storage
///
/// @param [in] pDef The Geant4 particle that will be added
void addParticle( G4ParticleDefinition* pDef);
    
    /// The internal storage consisting of PDG ID and the Geant4 particle
    std::unordered_map<int,G4ParticleDefinition*> m_pdgG4ParticleMap;
  };
}