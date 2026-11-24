// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Manager.hpp"

#include "ActsExamples/Geant4/MaterialPhysicsList.hpp"
#include "ActsExamples/Geant4/PhysicsListFactory.hpp"
#include "ActsExamples/Geant4/RunManager.hpp"

#include <memory>
#include <stdexcept>

#include <FTFP_BERT.hh>
#include <FTFP_BERT_ATL.hh>
#include <G4EmParameters.hh>
#include <G4FieldManager.hh>
#include <G4HadronicParameters.hh>
#include <G4HadronicProcessStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4MagneticField.hh>
#include <G4RunManager.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPhysicsList.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4Version.hh>

namespace ActsExamples {

Geant4Handle::Geant4Handle(std::unique_ptr<G4RunManager> _runManager,
                           std::unique_ptr<G4VUserPhysicsList> _physicsList,
                           std::string _physicsListName)
    : runManager(std::move(_runManager)),
      physicsList(_physicsList.get()),
      physicsListName(std::move(_physicsListName)) {
  if (runManager == nullptr) {
    throw std::invalid_argument("runManager cannot be null");
  }
  if (physicsList == nullptr) {
    throw std::invalid_argument("physicsList cannot be null");
  }

  // Set physics list
  runManager->SetUserInitialization(_physicsList.release());
}

Geant4Handle::~Geant4Handle() {
  clearMagneticField();
}

void Geant4Handle::clearMagneticField() {
  if (!m_fieldManager) {
    return;
  }
  // Multiple detector constructions may reuse logical volumes. Remove only
  // the field manager that this handle installed, wherever it is still used.
  for (auto* volume : *G4LogicalVolumeStore::GetInstance()) {
    if (volume->GetFieldManager() == m_fieldManager.get()) {
      volume->AssignFieldManager(nullptr);
    }
  }
  m_fieldManager.reset();
  m_magneticField.reset();
}

void Geant4Handle::tweakLogging(int level) const {
  Geant4Manager::tweakLogging(*runManager, level);
}

Geant4Manager& Geant4Manager::instance() {
  static Geant4Manager manager;
  return manager;
}

void Geant4Manager::tweakLogging(G4RunManager& runManager, int level) {
  runManager.SetVerboseLevel(level);
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

std::shared_ptr<Geant4Handle> Geant4Manager::currentHandle() const {
  std::lock_guard lock(m_mutex);
  return m_handle.lock();
}

std::shared_ptr<Geant4Handle> Geant4Manager::createHandle(
    const std::string& physicsList) {
  reserveRunManager();
  std::unique_ptr<G4VUserPhysicsList> list;
  try {
    list = createPhysicsList(physicsList);
  } catch (...) {
    // No run manager exists yet. A different factory can still be tried.
    std::lock_guard lock(m_mutex);
    m_created = false;
    throw;
  }
  return makeHandle(std::move(list), physicsList);
}

std::shared_ptr<Geant4Handle> Geant4Manager::createHandle(
    std::unique_ptr<G4VUserPhysicsList> physicsList,
    std::string physicsListName) {
  if (!physicsList) {
    throw std::invalid_argument("physicsList cannot be null");
  }
  reserveRunManager();
  return makeHandle(std::move(physicsList), std::move(physicsListName));
}

void Geant4Manager::reserveRunManager() {
  std::lock_guard lock(m_mutex);
  if (m_created) {
    throw std::runtime_error(
        "Geant4 run manager already created in this process");
  }
  m_created = true;
}

std::shared_ptr<Geant4Handle> Geant4Manager::makeHandle(
    std::unique_ptr<G4VUserPhysicsList> physicsList,
    std::string physicsListName) {
  auto runManager = std::make_unique<Geant4::RunManager>();

  auto handle = std::make_shared<Geant4Handle>(std::move(runManager),
                                               std::move(physicsList),
                                               std::move(physicsListName));

  std::lock_guard lock(m_mutex);
  m_handle = handle;
  return handle;
}

void Geant4Manager::registerPhysicsListFactory(
    std::string name,
    std::shared_ptr<Geant4::PhysicsListFactory> physicsListFactory) {
  if (!physicsListFactory) {
    throw std::invalid_argument("Physics list factory cannot be null");
  }
  std::lock_guard lock(m_mutex);
  if (m_physicsListFactories.contains(name)) {
    throw std::invalid_argument("name already mapped");
  }
  m_physicsListFactories.emplace(std::move(name),
                                 std::move(physicsListFactory));
}

std::unique_ptr<G4VUserPhysicsList> Geant4Manager::createPhysicsList(
    const std::string& name) const {
  std::shared_ptr<Geant4::PhysicsListFactory> factory;
  {
    std::lock_guard lock(m_mutex);
    const auto it = m_physicsListFactories.find(name);
    if (it == m_physicsListFactories.end()) {
      throw std::invalid_argument("name not mapped");
    }
    factory = it->second;
  }
  auto physics = factory->factorize();
  if (!physics) {
    throw std::invalid_argument("Physics list factory returned null");
  }
  return physics;
}

const std::unordered_map<std::string,
                         std::shared_ptr<Geant4::PhysicsListFactory>>&
Geant4Manager::getPhysicsListFactories() const {
  return m_physicsListFactories;
}

Geant4Manager::Geant4Manager() {
  registerPhysicsListFactory(
      "FTFP_BERT", std::make_shared<Geant4::PhysicsListFactoryFunction>(
                       []() { return std::make_unique<FTFP_BERT>(); }));
  registerPhysicsListFactory(
      "FTFP_BERT_ATL", std::make_shared<Geant4::PhysicsListFactoryFunction>(
                           []() { return std::make_unique<FTFP_BERT_ATL>(); }));
  registerPhysicsListFactory(
      "MaterialPhysicsList",
      std::make_shared<Geant4::PhysicsListFactoryFunction>(
          []() { return std::make_unique<Geant4::MaterialPhysicsList>(); }));
}

Geant4Manager::~Geant4Manager() = default;

}  // namespace ActsExamples
