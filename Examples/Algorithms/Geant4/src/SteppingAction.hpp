// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/MaterialInteractor.hpp"
#include "ActsExamples/EventData/SimHit.hpp"

#include <vector>

#include <G4UserSteppingAction.hh>
#include <globals.hh>

namespace ActsExamples::Geant4 {

/// @class SteppingAction
///
/// @brief Collects the RecordedMaterialSlab entities
///
/// The SteppingAction class is the implementation of the
/// Geant4 class SteppingAction. It extracts the weighted material
/// of every step and collects all material steps.
class SteppingAction final : public G4UserSteppingAction {
 public:
  /// Static access method to the instance
  static SteppingAction* instance();

  /// Construct the action and ensure singleton usage.
  SteppingAction();
  ~SteppingAction() final override;

  /// @brief Interface Method doing the step
  /// @note it creates and collects the MaterialInteraction entities
  /// @param step is the Geant4 step of the particle
  void UserSteppingAction(const G4Step* step) final override;

  /// Clear the recorded steps.
  void clear();

  /// Access the recorded material steps.
  const std::vector<Acts::MaterialInteraction>& materialSteps() const {
    return m_materialSteps;
  }

  /// Access the recorded tracking steps.
  const ActsExamples::SimHitContainer::sequence_type& trackSteps() const {
    return m_trackSteps;
  }

 private:
  /// Instance of the SteppingAction
  static SteppingAction* s_instance;

  /// The collected Acts::MaterialInteraction entities
  std::vector<Acts::MaterialInteraction> m_materialSteps;
  /// The collected ActsExamples::SimHit entities
  ActsExamples::SimHitContainer::sequence_type m_trackSteps;
};

}  // namespace ActsExamples::Geant4
