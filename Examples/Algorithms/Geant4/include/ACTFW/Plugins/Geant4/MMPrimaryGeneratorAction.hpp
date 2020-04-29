// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MMPrimaryGeneratorAction.hpp
///////////////////////////////////////////////////////////////////

#pragma once

#include <memory>

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

namespace FW {
namespace Geant4 {

/// @class MMPrimaryGeneratorAction
///
/// @brief configures the run
///
/// The MMPrimaryGeneratorAction is the implementation of the Geant4
/// class G4VUserPrimaryGeneratorAction. It generates a random direction
/// and shoots a geantino.
///
/// @todo tempate with RandomService
class MMPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  /// Constructor
  MMPrimaryGeneratorAction(const G4String& particleName = "geantino",
                           G4double energy = 1000. * MeV,
                           G4int randomSeed1 = 12345,
                           G4int randomSeed2 = 23456);
  /// Destructor
  ~MMPrimaryGeneratorAction() override;

  /// Static access method
  static MMPrimaryGeneratorAction* Instance();

  /// Interface method to generate the primary
  void GeneratePrimaries(G4Event*) final override;

  /// Access method to get the initial direction
  G4ThreeVector direction() { return m_direction; }

  /// Access method to get the start position
  G4ThreeVector startPosition() { return m_position; }

 private:
  /// Instance of the PrimaryGeneratorAction
  static MMPrimaryGeneratorAction* fgInstance;

  /// Pointer to the G4 particle gun
  std::unique_ptr<G4ParticleGun> fParticleGun;

  /// position to be returned
  G4ThreeVector m_position;
  /// direction to be returned
  G4ThreeVector m_direction;
};

}  // namespace Geant4
}  // namespace FW
