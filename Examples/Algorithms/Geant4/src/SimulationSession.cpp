// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "SimulationSession.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/detail/ScopeExit.hpp"
#include "ActsExamples/Geant4/MagneticFieldWrapper.hpp"
#include "ActsPlugins/FpeMonitoring/FpeMonitor.hpp"

#include <cmath>
#include <mutex>
#include <stdexcept>
#include <utility>

#include <G4FieldManager.hh>
#include <G4LogicalVolume.hh>
#include <G4PropagatorInField.hh>
#include <G4SystemOfUnits.hh>
#include <G4TransportationManager.hh>
#include <G4VPhysicalVolume.hh>
#include <Randomize.hh>

namespace ActsExamples::Geant4 {
namespace {
RunManager& checkedManager(const Geant4Handle& handle,
                           const std::string& physicsList) {
  if (handle.physicsListName != physicsList) {
    throw std::runtime_error("inconsistent physics list");
  }
  auto* manager = dynamic_cast<RunManager*>(handle.runManager.get());
  if (manager == nullptr) {
    throw std::invalid_argument(
        "Geant4 simulation requires ActsExamples::Geant4::RunManager");
  }
  manager->checkThread();
  if (manager->isRunActive()) {
    throw std::logic_error("Cannot reconfigure an active Geant4 run");
  }
  return *manager;
}

// The propagator is shared by successive configurations on this thread.
void processEventWithStepLimit(RunManager& manager, double largestStep) {
  if (!std::isfinite(largestStep)) {
    manager.processEvent();
    return;
  }
  auto* propagator = G4TransportationManager::GetTransportationManager()
                         ->GetPropagatorInField();
  const double previousStep = propagator->GetLargestAcceptableStep();
  Acts::detail::ScopeExit restore([propagator, previousStep] {
    propagator->SetLargestAcceptableStep(previousStep);
  });
  propagator->SetLargestAcceptableStep(largestStep / Acts::UnitConstants::mm *
                                       CLHEP::mm);
  manager.processEvent();
}
}  // namespace

SimulationSession::SimulationSession(std::shared_ptr<Geant4Handle> geant4Handle,
                                     const std::string& physicsList,
                                     double largestStep)
    : handle(std::move(geant4Handle)),
      m_manager(checkedManager(*handle, physicsList)),
      m_largestStep(largestStep) {}

void SimulationSession::install() {
  std::lock_guard lock(handle->mutex);
  // Configurations may share the handle sequentially. Geant4 retains only
  // the latest geometry and actions, so earlier sessions cannot start a run.
  if (m_manager.isRunActive()) {
    throw std::logic_error("Cannot reconfigure an active Geant4 run");
  }
  handle->clearMagneticField();
  handle->m_installedContext = eventContext;
}

void SimulationSession::setMagneticField(
    G4VPhysicalVolume& world,
    std::shared_ptr<const Acts::MagneticFieldProvider> field) {
  MagneticFieldWrapper::Config cfg;
  cfg.magneticField = std::move(field);
  handle->m_magneticField = std::make_unique<MagneticFieldWrapper>(cfg);
  handle->m_fieldManager = std::make_unique<G4FieldManager>();
  handle->m_fieldManager->SetDetectorField(handle->m_magneticField.get());
  handle->m_fieldManager->CreateChordFinder(handle->m_magneticField.get());
  world.GetLogicalVolume()->SetFieldManager(handle->m_fieldManager.get(), true);
}

void SimulationSession::beginRun() {
  std::lock_guard lock(handle->mutex);
  m_manager.checkThread();
  if (handle->m_installedContext != eventContext) {
    throw std::logic_error(
        "This Geant4 session was replaced on the shared handle");
  }
  ActsPlugins::FpeMonitor mon{0};
  if (configureActions) {
    if (m_manager.isRunActive()) {
      throw std::logic_error("Cannot configure actions during an active run");
    }
    configureActions();
    configureActions = {};
  }
  m_manager.beginRun();
  m_runStarted = true;
}

EventResult SimulationSession::simulate(const SimParticleContainer& particles,
                                        const Acts::GeometryContext& geometry,
                                        std::uint64_t seed) {
  std::lock_guard lock(handle->mutex);
  if (!m_runStarted) {
    throw std::logic_error("Geant4 run has not been initialized");
  }
  m_manager.checkThread();
  G4Random::setTheSeed(seed);
  EventStore event(geometry, particles);
  auto binding = eventContext->bind(event);
  ActsPlugins::FpeMonitor mon{0};
  try {
    processEventWithStepLimit(m_manager, m_largestStep);
  } catch (...) {
    m_runStarted = m_manager.isRunActive();
    throw;
  }
  return std::move(event.result);
}

void SimulationSession::endRun() {
  std::lock_guard lock(handle->mutex);
  if (m_runStarted) {
    m_manager.checkThread();
    m_runStarted = false;
    ActsPlugins::FpeMonitor mon{0};
    m_manager.endRun();
  }
}

}  // namespace ActsExamples::Geant4
