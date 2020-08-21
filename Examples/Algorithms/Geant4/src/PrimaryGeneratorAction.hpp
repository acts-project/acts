// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <globals.hh>

class G4ParticleGun;
class G4Event;

namespace ActsExamples {

/// @class PrimaryGeneratorAction
///
/// @brief configures the run
///
/// The PrimaryGeneratorAction is the implementation of the Geant4
/// class G4VUserPrimaryGeneratorAction. It generates a random direction
/// and shoots a geantino.
///
/// @todo tempate with RandomService
class PrimaryGeneratorAction final : public G4VUserPrimaryGeneratorAction {
 public:
  /// Static access method
  static PrimaryGeneratorAction* instance();

  /// Construct the action and ensure singleton usage.
  PrimaryGeneratorAction(const G4String& particleName = "geantino",
                         G4double energy = 1000. * MeV,
                         G4int randomSeed1 = 12345, G4int randomSeed2 = 23456);
  ~PrimaryGeneratorAction() final override;

  /// Interface method to generate the primary
  void GeneratePrimaries(G4Event*) final override;

  /// Access method to get the initial direction
  const G4ThreeVector& direction() const { return m_direction; }

  /// Access method to get the initial position
  const G4ThreeVector& position() const { return m_position; }

 private:
  /// Instance of the PrimaryGeneratorAction
  static PrimaryGeneratorAction* s_instance;

  /// Pointer to the G4 particle gun
  std::unique_ptr<G4ParticleGun> m_particleGun;
  /// position to be returned
  G4ThreeVector m_position;
  /// direction to be returned
  G4ThreeVector m_direction;
};

}  // namespace ActsExamples
