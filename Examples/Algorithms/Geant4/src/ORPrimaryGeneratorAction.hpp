// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

namespace ActsExamples {

  /// @class ORPrimaryGeneratorAction
  ///
  /// @brief configures the run
  ///
  /// The ORPrimaryGeneratorAction is the implementation of the Geant4
  /// class G4VUserPrimaryGeneratorAction. It generates a random direction
  /// and shoots a geantino.
  ///
  /// @todo tempate with RandomService
  class ORPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
  {
  public:
    /// Constructor
    ORPrimaryGeneratorAction(
                             G4int           randomSeed1  = 12345,
                             G4int           randomSeed2  = 23456);
    /// Destructor
    ~ORPrimaryGeneratorAction() override;

    /// Static access method
    static ORPrimaryGeneratorAction*
    instance();

    /// Interface method to generate the primary
    void
    GeneratePrimaries(G4Event*) final override;

    void prepareParticleGun(G4int pdg,
    G4double        momentum,
	 G4ThreeVector pos, G4ThreeVector dir);
	 
  private:
    /// Instance of the PrimaryGeneratorAction
    static ORPrimaryGeneratorAction* fgInstance;

    /// Pointer to the G4 particle gun
    std::unique_ptr<G4ParticleGun> fParticleGun;
  };

}  // namespace ActsExamples
