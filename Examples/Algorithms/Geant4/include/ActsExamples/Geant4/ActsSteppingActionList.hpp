// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"

#include <memory>
#include <optional>
#include <string>

#include <G4UserSteppingAction.hh>

namespace ActsExamples {

/// Geant4 only allows one user action of each type. This simple wrapper
/// dispatches multiple actions to Geant4.
class ActsSteppingActionList : public G4UserSteppingAction {
 public:
  struct Config {
    std::vector<G4UserSteppingAction *> actions;
  };

  ActsSteppingActionList(const Config &cfg) : m_cfg(cfg) {}

  void UserSteppingAction(const G4Step *step) override {
    for (auto action : m_cfg.actions) {
      if (action) {
        action->UserSteppingAction(step);
      }
    }
  }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
