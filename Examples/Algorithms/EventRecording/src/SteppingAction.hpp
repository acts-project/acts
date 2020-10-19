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

/// @class SteppingAction
///
/// @brief Collects the particles history
class SteppingAction : public G4UserSteppingAction {
 public:
  SteppingAction();
  ~SteppingAction() override;

  /// Static access method to the instance
  static SteppingAction* instance();

  /// @brief Interface Method doing the step and records the data
  /// @param step is the Geant4 step of the particle
  void UserSteppingAction(const G4Step* step) final override;

  /// Interface reset method
  void clear();

 private:
  /// Instance of the SteppingAction
  static SteppingAction* s_instance;
  /// The end vertex of the previous step
  std::shared_ptr<HepMC3::GenVertex> m_previousVertex = nullptr;
};
}  // namespace ActsExamples