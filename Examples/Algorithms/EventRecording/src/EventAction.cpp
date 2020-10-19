// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "EventAction.hpp"
#include <stdexcept>
#include <G4Event.hh>
#include <G4RunManager.hh>
#include "SteppingAction.hpp"

ActsExamples::EventAction* ActsExamples::EventAction::s_instance = nullptr;

ActsExamples::EventAction* ActsExamples::EventAction::instance() {
  // Static acces function via G4RunManager
  return s_instance;
}

ActsExamples::EventAction::EventAction() : G4UserEventAction() {
  if (s_instance) {
    throw std::logic_error("Attempted to duplicate a singleton");
  } else {
    s_instance = this;
  }
}

ActsExamples::EventAction::~EventAction() {
  s_instance = nullptr;
}

void ActsExamples::EventAction::BeginOfEventAction(const G4Event*) {
  SteppingAction::instance()->clear();
  m_event =
      std::make_shared<HepMC3::GenEvent>(HepMC3::Units::GEV, HepMC3::Units::MM);
}

void ActsExamples::EventAction::EndOfEventAction(const G4Event*) {}

void ActsExamples::EventAction::clear() {
  m_event = nullptr;
  SteppingAction::instance()->clear();
}

std::shared_ptr<HepMC3::GenEvent> ActsExamples::EventAction::event() const {
  return m_event;
}