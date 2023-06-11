// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Manager.hpp"

#include "ActsExamples/Geant4/PhysicsListFactory.hpp"

#include <memory>
#include <stdexcept>

#include <FTFP_BERT.hh>
#include <FTFP_BERT_ATL.hh>
#include <G4EmParameters.hh>
#include <G4HadronicParameters.hh>
#include <G4HadronicProcessStore.hh>
#include <G4Profiler.hh>
#include <G4RunManager.hh>
#include <G4RunManagerFactory.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPhysicsList.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4Version.hh>

namespace ActsExamples {

Geant4Instance::Geant4Instance(std::unique_ptr<G4RunManager> _runManager)
    : runManager(std::move(_runManager)) {}

Geant4Instance::~Geant4Instance() {
  clearPhysicsList();
}

void Geant4Instance::registerPhysicsList(
    std::string name, std::unique_ptr<G4VUserPhysicsList> physicsList) {
  if (physicsLists.find(name) != physicsLists.end()) {
    throw std::invalid_argument("name already mapped");
  }
  physicsLists.emplace(std::move(name), physicsList.release());
}

G4VUserPhysicsList* Geant4Instance::getPhysicsList(
    const std::string& name) const {
  auto it = physicsLists.find(name);
  if (it == physicsLists.end()) {
    throw std::invalid_argument("name not mapped");
  }
  return it->second;
}

G4VUserPhysicsList* Geant4Instance::createRegisterAndGetPhysicsList(
    const std::string& name) {
  auto it = physicsLists.find(name);
  if (it == physicsLists.end()) {
    auto physicsList = Geant4Manager::instance().createPhysicsList(name);
    it = physicsLists.emplace(name, physicsList.release()).first;
  }
  return it->second;
}

G4VUserPhysicsList* Geant4Instance::registerAndGetAnonymousPhysicsList(
    std::unique_ptr<G4VUserPhysicsList> physicsList) {
  delete anonymousPhysicList;
  anonymousPhysicList = physicsList.release();
  return anonymousPhysicList;
}

void Geant4Instance::tweekLogging(int level) const {
  runManager->SetVerboseLevel(level);
  G4EventManager::GetEventManager()->SetVerboseLevel(level);
  G4EventManager::GetEventManager()->GetTrackingManager()->SetVerboseLevel(
      level);
  G4EventManager::GetEventManager()->GetStackManager()->SetVerboseLevel(level);

  // Suppress the printing of physics information.
#if G4VERSION_NUMBER >= 1100
  G4HadronicParameters::Instance()->SetVerboseLevel(0);
  G4HadronicProcessStore::Instance()->SetVerbose(0);
  G4EmParameters::Instance()->SetIsPrintedFlag(true);
#endif
}

void Geant4Instance::clearPhysicsList() {
  // delete all physics lists except the active one since Geant4 will do that
  for (auto it = physicsLists.begin(); it != physicsLists.end();) {
    if (it->second != runManager->GetUserPhysicsList()) {
      delete it->second;
      it = physicsLists.erase(it);
    } else {
      ++it;
    }
  }
  if (anonymousPhysicList != runManager->GetUserPhysicsList()) {
    delete anonymousPhysicList;
  }
}

Geant4Manager& Geant4Manager::instance() {
  static Geant4Manager manager;
  return manager;
}

std::shared_ptr<Geant4Instance> Geant4Manager::create() {
  if (!m_instance.expired()) {
    throw std::runtime_error("creating a second instance is prohibited");
  }
  if (m_created) {
    throw std::runtime_error(
        "creating a new instance is prohibited. you have to hold onto the "
        "first one.");
  }

  auto runManager = std::unique_ptr<G4RunManager>(
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly));

  auto instance = std::make_shared<Geant4Instance>(std::move(runManager));

  m_created = true;
  m_instance = instance;
  return instance;
}

void Geant4Manager::registerPhysicsListFactory(
    std::string name, std::shared_ptr<PhysicsListFactory> physicsListFactory) {
  if (m_physicsListsFactory.find(name) != m_physicsListsFactory.end()) {
    throw std::invalid_argument("name already mapped");
  }
  m_physicsListsFactory.emplace(std::move(name), std::move(physicsListFactory));
}

std::unique_ptr<G4VUserPhysicsList> Geant4Manager::createPhysicsList(
    const std::string& name) const {
  auto it = m_physicsListsFactory.find(name);
  if (it == m_physicsListsFactory.end()) {
    throw std::invalid_argument("name not mapped");
  }
  return it->second->factorize();
}

Geant4Manager::Geant4Manager() {
  registerPhysicsListFactory(
      "FTFP_BERT", std::make_shared<PhysicsListFactoryFunction>(
                       []() { return std::make_unique<FTFP_BERT>(); }));
  registerPhysicsListFactory(
      "FTFP_BERT_ATL", std::make_shared<PhysicsListFactoryFunction>(
                           []() { return std::make_unique<FTFP_BERT_ATL>(); }));
}

Geant4Manager::~Geant4Manager() = default;

}  // namespace ActsExamples
