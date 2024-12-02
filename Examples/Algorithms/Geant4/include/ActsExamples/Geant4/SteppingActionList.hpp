// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <vector>

#include <G4UserSteppingAction.hh>

namespace ActsExamples::Geant4 {

/// Geant4 only allows one user action of each type. This simple wrapper
/// dispatches multiple actions to Geant4.
class SteppingActionList : public G4UserSteppingAction {
 public:
  struct Config {
    std::vector<std::shared_ptr<G4UserSteppingAction>> actions;
  };

  explicit SteppingActionList(const Config &cfg) : m_cfg(cfg) {}

  void UserSteppingAction(const G4Step *step) override {
    for (const auto &action : m_cfg.actions) {
      if (action) {
        action->UserSteppingAction(step);
      }
    }
  }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples::Geant4
