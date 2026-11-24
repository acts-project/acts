// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "NativeWorkerTransport.hpp"

#include <exception>
#include <memory>
#include <stdexcept>
#include <utility>

#include <G4Box.hh>
#include <G4Electron.hh>
#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4Gamma.hh>
#include <G4Geantino.hh>
#include <G4LogicalVolume.hh>
#include <G4MTRunManager.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4ParticleGun.hh>
#include <G4Positron.hh>
#include <G4Proton.hh>
#include <G4RunManager.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4UserEventAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4VUserActionInitialization.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserEventInformation.hh>
#include <G4VUserPhysicsList.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <Randomize.hh>

#include "NativeWorkerExchange.hpp"

namespace ActsTests::NativeWorkers {
namespace {

class Telescope final : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* world = new G4LogicalVolume(
        new G4Box("World", 100 * CLHEP::mm, 100 * CLHEP::mm, 100 * CLHEP::mm),
        nist->FindOrBuildMaterial("G4_Galactic"), "World");
    auto* placement =
        new G4PVPlacement(nullptr, {}, world, "World", nullptr, false, 0);
    auto* sensor = new G4LogicalVolume(
        new G4Box("Sensor", 90 * CLHEP::mm, 90 * CLHEP::mm, 0.25 * CLHEP::mm),
        nist->FindOrBuildMaterial("G4_Si"), "Sensor");
    new G4PVPlacement(nullptr, {0, 0, 30 * CLHEP::mm}, sensor, "Sensor", world,
                      false, 0);
    new G4PVPlacement(nullptr, {0, 0, 60 * CLHEP::mm}, sensor, "Sensor", world,
                      false, 1);
    return placement;
  }
};

class Physics final : public G4VUserPhysicsList {
 public:
  void ConstructParticle() override {
    G4Geantino::GeantinoDefinition();
    // Geant4's default production cuts reference these four definitions.
    G4Gamma::Definition();
    G4Electron::Definition();
    G4Positron::Definition();
    G4Proton::Definition();
  }
  void ConstructProcess() override { AddTransportation(); }
  void SetCuts() override { SetCutsWithDefault(); }
};

// Geant4 owns this object for exactly one event on one worker. Keep the ACTS
// request identity separate from Geant4's native event scheduling counter.
struct EventData final : G4VUserEventInformation {
  explicit EventData(std::shared_ptr<Job> input)
      : job(std::move(input)), result{job->request.eventId, 0., {}} {}

  void Print() const override {}
  std::shared_ptr<Job> job;
  Result result;
  std::exception_ptr error;
};

class Generator final : public G4VUserPrimaryGeneratorAction {
 public:
  explicit Generator(Exchange& exchange) : m_exchange(exchange) {
    m_gun.SetParticleDefinition(G4Geantino::GeantinoDefinition());
    m_gun.SetParticleMomentumDirection({0, 0, 1});
    m_gun.SetParticleEnergy(CLHEP::GeV);
  }

  void GeneratePrimaries(G4Event* event) override {
    auto job = m_exchange.take();
    if (!job) {
      // End this worker's run after its current, empty event. Every worker
      // observes closure independently, including workers idle in take().
      G4RunManager::GetRunManager()->AbortRun(true);
      return;
    }
    auto* data = new EventData(std::move(job));
    event->SetUserInformation(data);
    try {
      G4Random::setTheSeed(data->job->request.seed);
      data->result.random = G4UniformRand();
      const double x = data->job->request.x + data->result.random;
      m_gun.SetParticlePosition({x * CLHEP::mm, 0, 0});
      m_gun.GeneratePrimaryVertex(event);
    } catch (...) {
      data->error = std::current_exception();
    }
  }

 private:
  Exchange& m_exchange;
  G4ParticleGun m_gun{1};
};

class Stepping final : public G4UserSteppingAction {
 public:
  explicit Stepping(std::function<void(int)> callback)
      : m_onTransportStep(std::move(callback)) {}

  void UserSteppingAction(const G4Step* step) override {
    const auto* pre = step->GetPreStepPoint();
    if (pre->GetStepStatus() != fGeomBoundary ||
        pre->GetPhysicalVolume()->GetName() != "Sensor") {
      return;
    }
    const auto* event =
        G4EventManager::GetEventManager()->GetConstCurrentEvent();
    auto* data = static_cast<EventData*>(event->GetUserInformation());
    try {
      const auto& position = pre->GetPosition();
      data->result.hits.emplace_back(position.x() / CLHEP::mm,
                                     position.z() / CLHEP::mm);
      if (data->result.hits.size() == 1) {
        m_onTransportStep(data->job->request.eventId);
        if (data->job->request.fail) {
          throw TransportFailure("Injected Geant4 transport failure");
        }
      }
    } catch (...) {
      data->error = std::current_exception();
      G4RunManager::GetRunManager()->AbortEvent();
    }
  }

 private:
  std::function<void(int)> m_onTransportStep;
};

class Publish final : public G4UserEventAction {
 public:
  explicit Publish(Exchange& exchange) : m_exchange(exchange) {}
  void EndOfEventAction(const G4Event* event) override {
    auto* data = static_cast<EventData*>(event->GetUserInformation());
    if (data == nullptr) {
      return;  // Empty terminal slot.
    }
    if (!data->error && event->IsAborted()) {
      data->error =
          std::make_exception_ptr(std::runtime_error("Geant4 event aborted"));
    }
    if (data->error) {
      m_exchange.fail(data->job, data->error);
    } else {
      m_exchange.complete(data->job, std::move(data->result));
    }
  }

 private:
  Exchange& m_exchange;
};

class Actions final : public G4VUserActionInitialization {
 public:
  Actions(Exchange& exchange, std::function<void(int)> callback)
      : m_exchange(exchange), m_onTransportStep(std::move(callback)) {}
  void Build() const override {
    SetUserAction(new Generator(m_exchange));
    SetUserAction(new Stepping(m_onTransportStep));
    SetUserAction(new Publish(m_exchange));
  }

 private:
  Exchange& m_exchange;
  std::function<void(int)> m_onTransportStep;
};

}  // namespace

void configure(G4MTRunManager& manager, Exchange& exchange,
               std::function<void(int)> onTransportStep) {
  manager.SetEventModulo(1);
  manager.SetUserInitialization(new Telescope());
  manager.SetUserInitialization(new Physics());
  manager.SetUserInitialization(
      new Actions(exchange, std::move(onTransportStep)));
}

}  // namespace ActsTests::NativeWorkers
