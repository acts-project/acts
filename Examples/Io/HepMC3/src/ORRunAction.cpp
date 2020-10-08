// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ORRunAction.hpp"

#include <stdexcept>

#include <G4Run.hh>

#include "OREventAction.hpp"

using namespace ActsExamples;

ORRunAction* ORRunAction::s_instance = nullptr;

ORRunAction* ORRunAction::instance() {
  return s_instance;
}

ORRunAction::ORRunAction() : G4UserRunAction() {
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate the ORRunAction singleton");
  } else {
    s_instance = this;
  }
}

ORRunAction::~ORRunAction() {
  s_instance = nullptr;
}

void ORRunAction::BeginOfRunAction(const G4Run* aRun) {
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  // initialize event cumulative quantities
  OREventAction::instance()->clear();
}

void ORRunAction::EndOfRunAction(const G4Run* aRun) {
  G4int nofEvents = aRun->GetNumberOfEvent();
  if (nofEvents == 0)
    return;

  // Print
  G4cout << "\n--------------------End of Run------------------------------\n"
         << "\n------------------------------------------------------------\n"
         << G4endl;
}
