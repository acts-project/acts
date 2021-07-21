// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "RunAction.hpp"

#include <stdexcept>

#include <G4Run.hh>

#include "EventAction.hpp"

namespace ActsExamples::Geant4 {

RunAction* RunAction::s_instance = nullptr;

RunAction* RunAction::instance() {
  return s_instance;
}

RunAction::RunAction() : G4UserRunAction() {
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate the RunAction singleton");
  } else {
    s_instance = this;
  }
}

RunAction::~RunAction() {
  s_instance = nullptr;
}

void RunAction::BeginOfRunAction(const G4Run* aRun) {
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  // initialize event cumulative quantities
  EventAction::instance()->clear();
}

void RunAction::EndOfRunAction(const G4Run* aRun) {
  G4int nofEvents = aRun->GetNumberOfEvent();
  if (nofEvents == 0)
    return;

  // Print
  G4cout << "\n--------------------End of Run------------------------------\n"
         << "\n------------------------------------------------------------\n"
         << G4endl;
}

}  // namespace ActsExamples::Geant4