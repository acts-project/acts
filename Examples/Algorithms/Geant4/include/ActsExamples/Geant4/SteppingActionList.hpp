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

  /// Register a batch of actions on the installed list before a run starts.
  /// @param extensions New actions in callback order.
  void appendActions(const Config &extensions) {
    for (const auto &action : extensions.actions) {
      if (action != nullptr) {
        action->SetSteppingManagerPointer(fpSteppingManager);
      }
    }
    m_cfg.actions.insert(m_cfg.actions.end(), extensions.actions.begin(),
                         extensions.actions.end());
  }

  /// Forward the manager to the delegated actions when Geant4 installs this
  /// list.
  /// @param manager Geant4 stepping manager for the current worker.
  void SetSteppingManagerPointer(G4SteppingManager *manager) override {
    G4UserSteppingAction::SetSteppingManagerPointer(manager);
    for (const auto &action : m_cfg.actions) {
      if (action != nullptr) {
        action->SetSteppingManagerPointer(manager);
      }
    }
  }

  void UserSteppingAction(const G4Step *stepPtr) override {
    for (const auto &action : m_cfg.actions) {
      if (action != nullptr) {
        action->UserSteppingAction(stepPtr);
      }
    }
  }

 private:
  Config m_cfg;
};

}  // namespace ActsExamples::Geant4
