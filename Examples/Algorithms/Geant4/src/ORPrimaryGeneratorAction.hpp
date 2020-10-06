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
    //~ ORPrimaryGeneratorAction(const G4int& pdg = 211,
                             //~ G4double        momentum       = 1000. * MeV,
                             //~ G4bool lockAngle = false,
                             //~ G4double phi = 0.,
                             //~ G4double theta = 0.5 * M_PI,
                             //~ G4bool lockPosition = true,
                             //~ G4ThreeVector pos = {0., 0., 0.},
                             //~ G4int           randomSeed1  = 12345,
                             //~ G4int           randomSeed2  = 23456);
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

    //~ /// Access method to get the initial direction
    //~ G4ThreeVector
    //~ direction()
    //~ {
      //~ return m_direction;
    //~ }

    //~ /// Access method to get the start position
    //~ G4ThreeVector
    //~ startPosition()
    //~ {
      //~ return m_position;
    //~ }

    void prepareParticleGun(G4int pdg,
    G4double        momentum,
	 G4double phi,
	 G4double theta,
	 G4ThreeVector pos);
	 
  private:
    /// Instance of the PrimaryGeneratorAction
    static ORPrimaryGeneratorAction* fgInstance;

    /// Pointer to the G4 particle gun
    std::unique_ptr<G4ParticleGun> fParticleGun;

    //~ /// position to be returned
    //~ G4ThreeVector m_position;
    //~ /// direction to be returned
    //~ G4ThreeVector m_direction;
    //~ G4bool m_lockAngle = false;
	//~ G4double m_phi = 0.;
	//~ G4double m_theta = 0.5 * M_PI;
	//~ G4bool m_lockPosition = true;
	//~ G4ThreeVector m_pos = {0., 0., 0.};
  };

}  // namespace ActsExamples
