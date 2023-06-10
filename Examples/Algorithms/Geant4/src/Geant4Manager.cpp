// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Manager.hpp"

#include <memory>
#include <stdexcept>

#include <G4RunManager.hh>
#include <G4RunManagerFactory.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPhysicsList.hh>
#include <G4VUserPrimaryGeneratorAction.hh>

namespace ActsExamples {

Geant4Instance::Geant4Instance(G4RunManager* _runManager)
    : runManager(_runManager) {}

Geant4Instance::~Geant4Instance() = default;

Geant4Manager& Geant4Manager::instance() {
  static Geant4Manager manager;
  return manager;
}

std::shared_ptr<Geant4Instance> Geant4Manager::create() {
  if (!m_instance.expired()) {
    throw std::runtime_error("creating a second instance is prohibited");
  }

  auto instance = std::make_shared<Geant4Instance>(m_runManager);
  m_instance = instance;
  return instance;
}

Geant4Manager::Geant4Manager() {
  // proper memory management might cause segfaults
  m_runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
}

Geant4Manager::~Geant4Manager() = default;

}  // namespace ActsExamples
