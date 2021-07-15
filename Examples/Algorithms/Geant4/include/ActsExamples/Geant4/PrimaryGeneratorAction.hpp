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

namespace ActsExamples::Geant4 {

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
  struct Config {
    /// Name of the generated particle
    G4String particleName = "geantino";
    /// Energy of the generated particles
    G4double energy = 1000.0;
    /// First random seed
    G4int randomSeed1 = 12345;
    /// Second random seed
    G4int randomSeed2 = 45678;
    /// Position in X of the generated vertex
    G4double vertexPosX = 0.0;
    /// Position in Y of the generated vertex
    G4double vertexPosY = 0.0;
    /// Position in Z of the generated vertex
    G4double vertexPosZ = 0.0;
    /// Spread in X of the generated vertex
    G4double vertexSigmaX = 0.0;
    /// Spread in Y of the generated vertex
    G4double vertexSigmaY = 0.0;
    /// Spread in Z of the generated vertex
    G4double vertexSigmaZ = 0.0;
    /// Azimutal angle phi range of the generated particles
    std::pair<G4double, G4double> phiRange = {-M_PI, M_PI};
    /// Pseudorapidity eta range of the generated particles
    std::pair<G4double, G4double> etaRange = {-4., 4.};
    /// Variable from which the particle generation is uniform.
    G4String samplingVariable = "theta";
  };

  /// Static access method
  static PrimaryGeneratorAction* instance();

  /// Construct the action and ensure singleton usage.
  PrimaryGeneratorAction(const Config& cfg);
  ~PrimaryGeneratorAction() final override;

  /// Interface method to generate the primary
  void GeneratePrimaries(G4Event*) final override;

  /// Access method to get the initial direction
  const G4ThreeVector& direction() const { return m_direction; }

  /// Access method to get the initial position
  const G4ThreeVector& position() const { return m_position; }

 private:
  /// The config class
  Config m_cfg;

  /// Instance of the PrimaryGeneratorAction
  static PrimaryGeneratorAction* s_instance;

  /// Pointer to the G4 particle gun
  std::unique_ptr<G4ParticleGun> m_particleGun;
  /// position to be returned
  G4ThreeVector m_position;
  /// direction to be returned
  G4ThreeVector m_direction;
};

}  // namespace ActsExamples::Geant4
