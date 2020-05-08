// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <G4UserSteppingAction.hh>
#include <globals.hh>
#include <vector>

#include "ACTFW/EventData/SimHit.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"

namespace ActsExamples {

/// @class SteppingAction
///
/// @brief Collects the RecordedMaterialProperties entities
///
/// The SteppingAction class is the implementation of the
/// Geant4 class SteppingAction. It extracts the weighted material
/// of every step and collects all material steps.
class SteppingAction final : public G4UserSteppingAction {
 public:
  /// Constructor
  SteppingAction();

  /// Destructor
  ~SteppingAction() final override;

  /// Static access method to the instance
  static SteppingAction* Instance();

  /// @brief Interface Method doing the step
  /// @note it creates and collects the MaterialInteraction entities
  /// @param step is the Geant4 step of the particle
  void UserSteppingAction(const G4Step* step) final override;

  /// Interface reset method
  /// @note it clears the collected step vector
  void Reset();

  /// Access to the collected Acts::MaterialInteraction entities
  std::vector<Acts::MaterialInteraction> materialSteps() { return m_steps; }

  /// Access to the collected FW::SimHitContainer entities
  FW::SimHitContainer::sequence_type trackSteps() { return m_tracksteps; }

 private:
  /// Instance of the SteppingAction
  static SteppingAction* s_instance;

  /// The collected Acts::MaterialInteraction entities
  std::vector<Acts::MaterialInteraction> m_steps = {};
  /// The collected FW::SimHit entities
  FW::SimHitContainer::sequence_type m_tracksteps;
};

}  // namespace ActsExamples
