// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/Geant4/RunManager.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <G4Box.hh>
#include <G4Electron.hh>
#include <G4Event.hh>
#include <G4Gamma.hh>
#include <G4Geantino.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4MaterialCutsCouple.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4ParticleGun.hh>
#include <G4Positron.hh>
#include <G4ProductionCuts.hh>
#include <G4ProductionCutsTable.hh>
#include <G4Proton.hh>
#include <G4Region.hh>
#include <G4Run.hh>
#include <G4StateManager.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPhysicsList.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <Randomize.hh>

namespace {

class CountingRunManager final : public ActsExamples::Geant4::RunManager {
 public:
  unsigned int initializations = 0;
  unsigned int terminations = 0;
  unsigned int deletions = 0;

  void RunInitialization() override {
    ++initializations;
    G4RunManager::RunInitialization();
  }
  void RunTermination() override {
    ++terminations;
    G4RunManager::RunTermination();
  }

 protected:
  void ReportEventDeletion(const G4Event* event) override {
    ++deletions;
    G4RunManager::ReportEventDeletion(event);
  }
};

class TestPhysicsList final : public G4VUserPhysicsList {
 public:
  void ConstructParticle() override {
    G4Geantino::GeantinoDefinition();
    G4Gamma::GammaDefinition();
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4Proton::ProtonDefinition();
  }
  void ConstructProcess() override { AddTransportation(); }
  void SetCuts() override { SetCutsWithDefault(); }
};

class Detector final : public G4VUserDetectorConstruction {
 public:
  explicit Detector(G4ProductionCuts& cuts) : m_cuts(cuts) {}

  G4VPhysicalVolume* Construct() override {
    if (world != nullptr) {
      return world;
    }
    auto* nist = G4NistManager::Instance();
    auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    silicon = nist->FindOrBuildMaterial("G4_Si");
    aluminium = nist->FindOrBuildMaterial("G4_Al");
    auto* worldLV = new G4LogicalVolume(
        new G4Box("world", CLHEP::m, CLHEP::m, CLHEP::m), vacuum, "world");
    world = new G4PVPlacement(nullptr, {}, worldLV, "world", nullptr, false, 0);
    slab = new G4LogicalVolume(
        new G4Box("slab", 10 * CLHEP::mm, 100 * CLHEP::mm, 100 * CLHEP::mm),
        silicon, "slab");
    new G4PVPlacement(nullptr, {}, slab, "slab", worldLV, false, 0);
    region = new G4Region("slab");
    region->SetProductionCuts(&m_cuts);
    region->AddRootLogicalVolume(slab);
    return world;
  }

  G4VPhysicalVolume* world = nullptr;
  G4LogicalVolume* slab = nullptr;
  G4Region* region = nullptr;
  G4Material* silicon = nullptr;
  G4Material* aluminium = nullptr;

 private:
  G4ProductionCuts& m_cuts;
};

struct Observations {
  std::vector<std::string> callbacks;
  std::vector<double> random;
  std::vector<std::pair<std::string, double>> steps;
  std::vector<int> eventIds;
  unsigned int begins = 0;
  unsigned int ends = 0;
  int planned = 0;
  int processed = 0;
  bool throwAtBegin = false;
  bool throwAtEnd = false;
  bool abortEvent = false;
  bool abortRun = false;
  bool softAbort = false;
};

class RunAction final : public G4UserRunAction {
 public:
  explicit RunAction(Observations& data) : m_data(data) {}
  G4Run* GenerateRun() override {
    m_data.callbacks.emplace_back("generate run");
    return new G4Run();
  }
  void BeginOfRunAction(const G4Run* run) override {
    BOOST_CHECK_EQUAL(G4StateManager::GetStateManager()->GetCurrentState(),
                      G4State_GeomClosed);
    m_data.callbacks.emplace_back("begin run");
    ++m_data.begins;
    m_data.planned = run->GetNumberOfEventToBeProcessed();
  }
  void EndOfRunAction(const G4Run* run) override {
    BOOST_CHECK_EQUAL(G4StateManager::GetStateManager()->GetCurrentState(),
                      G4State_GeomClosed);
    m_data.callbacks.emplace_back("end run");
    ++m_data.ends;
    m_data.planned = run->GetNumberOfEventToBeProcessed();
    m_data.processed = run->GetNumberOfEvent();
    if (m_data.throwAtEnd) {
      throw std::runtime_error("end run failure");
    }
  }

 private:
  Observations& m_data;
};

class Generator final : public G4VUserPrimaryGeneratorAction {
 public:
  explicit Generator(Observations& data) : m_data(data) {}
  void GeneratePrimaries(G4Event* event) override {
    m_data.callbacks.emplace_back("generate event");
    const double sample = G4UniformRand();
    m_data.random.push_back(sample);
    G4ParticleGun gun(G4Geantino::GeantinoDefinition());
    const double phi = 0.01 * sample;
    gun.SetParticleMomentumDirection({std::cos(phi), std::sin(phi), 0});
    gun.SetParticlePosition({-100 * CLHEP::mm, 0, 0});
    gun.SetParticleEnergy(CLHEP::GeV);
    gun.GeneratePrimaryVertex(event);
  }

 private:
  Observations& m_data;
};

class EventAction final : public G4UserEventAction {
 public:
  explicit EventAction(Observations& data) : m_data(data) {}
  void BeginOfEventAction(const G4Event* event) override {
    m_data.callbacks.emplace_back("begin event");
    m_data.random.push_back(G4UniformRand());
    m_data.eventIds.push_back(event->GetEventID());
    if (m_data.throwAtBegin) {
      throw std::runtime_error("event failure");
    }
    if (m_data.abortEvent) {
      G4RunManager::GetRunManager()->AbortEvent();
    }
    if (m_data.abortRun) {
      G4RunManager::GetRunManager()->AbortRun(m_data.softAbort);
    }
  }
  void EndOfEventAction(const G4Event*) override {
    m_data.callbacks.emplace_back("end event");
    m_data.random.push_back(G4UniformRand());
  }

 private:
  Observations& m_data;
};

class SteppingAction final : public G4UserSteppingAction {
 public:
  explicit SteppingAction(Observations& data) : m_data(data) {}
  void UserSteppingAction(const G4Step* step) override {
    m_data.steps.emplace_back(step->GetPreStepPoint()->GetMaterial()->GetName(),
                              step->GetStepLength());
  }

 private:
  Observations& m_data;
};

}  // namespace

BOOST_AUTO_TEST_SUITE(Geant4RunManager)

// A single lifecycle: Geant4 cannot safely construct a new manager per case.
BOOST_AUTO_TEST_CASE(PersistentRunLifecycle) {
  Observations data;
  G4ProductionCuts cuts;
  cuts.SetProductionCut(CLHEP::mm);
  G4ProductionCuts replacementCuts;
  replacementCuts.SetProductionCut(2 * CLHEP::mm);
  auto manager = std::make_unique<CountingRunManager>();
  auto* detector = new Detector(cuts);
  manager->SetUserInitialization(detector);
  manager->SetUserInitialization(new TestPhysicsList());
  manager->SetUserAction(new RunAction(data));
  manager->SetUserAction(new Generator(data));
  manager->SetUserAction(new EventAction(data));
  manager->SetUserAction(new SteppingAction(data));
  manager->Initialize();

  constexpr std::array<long, 3> seeds{17, 43, 101};
  std::vector<Observations> reference;
  for (long seed : seeds) {
    data = {};
    G4Random::setTheSeed(seed);
    manager->BeamOn(1);
    reference.push_back(data);
  }
  data = {};
  const auto initializations = manager->initializations;
  const auto deletions = manager->deletions;
  manager->beginRun();
  BOOST_CHECK_EQUAL(manager->initializations, initializations + 1);
  BOOST_CHECK_EQUAL(data.begins, 1);
  BOOST_CHECK_EQUAL(data.planned, 0);  // Unknown until ACTS supplies events.
  BOOST_CHECK_THROW(manager->beginRun(), std::logic_error);
  const auto* run = manager->GetCurrentRun();
  const auto* couple = detector->slab->GetMaterialCutsCouple();
  const std::vector<std::string> eventCallbacks = {"generate event",
                                                   "begin event", "end event"};
  for (std::size_t i = 0; i < seeds.size(); ++i) {
    data.callbacks.clear();
    data.random.clear();
    data.steps.clear();
    G4Random::setTheSeed(seeds[i]);
    manager->processEvent();
    BOOST_CHECK(data.callbacks == eventCallbacks);
    BOOST_CHECK(data.random == reference[i].random);
    BOOST_CHECK(data.steps == reference[i].steps);
    BOOST_CHECK_EQUAL(data.eventIds.back(), i);
    BOOST_CHECK_EQUAL(manager->GetCurrentRun(), run);
    BOOST_CHECK_EQUAL(detector->slab->GetMaterialCutsCouple(), couple);
    BOOST_CHECK_EQUAL(manager->initializations, initializations + 1);
    BOOST_CHECK_EQUAL(manager->deletions, deletions + i + 1);
    BOOST_CHECK_EQUAL(data.ends, 0);
  }
  manager->endRun();
  BOOST_CHECK_EQUAL(data.ends, 1);
  BOOST_CHECK_EQUAL(data.processed, 3);
  BOOST_CHECK_EQUAL(data.planned, 3);
  BOOST_CHECK_EQUAL(manager->GetCurrentRun()->GetRunID(), 3);
  BOOST_CHECK_EQUAL(G4StateManager::GetStateManager()->GetCurrentState(),
                    G4State_Idle);
  manager->endRun();
  BOOST_CHECK_EQUAL(data.ends, 1);
  BOOST_CHECK_THROW(manager->processEvent(), std::logic_error);
  BOOST_CHECK_THROW(manager->beginRun(-1), std::invalid_argument);

  BOOST_TEST_CONTEXT("zero events and explicit event count") {
    data = {};
    manager->beginRun();
    manager->endRun();
    BOOST_CHECK_EQUAL(data.begins, 1);
    BOOST_CHECK_EQUAL(data.ends, 1);
    BOOST_CHECK_EQUAL(data.processed, 0);
    BOOST_CHECK_EQUAL(data.planned, 0);
    data = {};
    manager->beginRun(2);
    BOOST_CHECK_EQUAL(data.planned, 2);
    manager->processEvent();
    manager->processEvent();
    manager->endRun();
    BOOST_CHECK_EQUAL(data.processed, 2);
    BOOST_CHECK_EQUAL(data.planned, 2);
  }

  BOOST_TEST_CONTEXT("reconfiguration between runs, including energy vectors") {
    detector->slab->SetMaterial(detector->aluminium);
    detector->region->SetProductionCuts(&replacementCuts);
    manager->GeometryHasBeenModified();
    manager->beginRun();
    BOOST_CHECK_EQUAL(detector->slab->GetMaterialCutsCouple()->GetMaterial(),
                      detector->aluminium);
    BOOST_CHECK_EQUAL(
        detector->slab->GetMaterialCutsCouple()->GetProductionCuts(),
        &replacementCuts);
    manager->processEvent();
    manager->endRun();
    auto* table = G4ProductionCutsTable::GetProductionCutsTable();
    table->SetEnergyCutVector(
        std::vector<double>(table->GetTableSize(), 0.7 * CLHEP::MeV), 0);
    manager->beginRun();
    for (double value : *table->GetEnergyCutsVector(0)) {
      BOOST_CHECK_EQUAL(value, 0.7 * CLHEP::MeV);
    }
    manager->processEvent();
    manager->endRun();
  }

  BOOST_TEST_CONTEXT("event abort permits the next event") {
    data = {};
    manager->beginRun();
    data.abortEvent = true;
    manager->processEvent();
    BOOST_CHECK(manager->isRunActive());
    data.abortEvent = false;
    manager->processEvent();
    manager->endRun();
    BOOST_CHECK_EQUAL(data.processed, 2);
    BOOST_CHECK_EQUAL(data.ends, 1);
  }

  for (bool soft : {false, true}) {
    BOOST_TEST_CONTEXT("run abort, soft=" << soft) {
      data = {};
      manager->beginRun();
      data.abortRun = true;
      data.softAbort = soft;
      BOOST_CHECK_THROW(manager->processEvent(), std::runtime_error);
      BOOST_CHECK(!manager->isRunActive());
      BOOST_CHECK_EQUAL(data.ends, 1);
      BOOST_CHECK_EQUAL(G4StateManager::GetStateManager()->GetCurrentState(),
                        G4State_Idle);
    }
  }

  BOOST_TEST_CONTEXT("event exception releases the current event and run") {
    data = {};
    manager->beginRun();
    const auto before = manager->deletions;
    data.throwAtBegin = true;
    BOOST_CHECK_THROW(manager->processEvent(), std::runtime_error);
    BOOST_CHECK(!manager->isRunActive());
    BOOST_CHECK_EQUAL(manager->GetCurrentEvent(), nullptr);
    BOOST_CHECK_EQUAL(manager->deletions, before + 1);
    BOOST_CHECK_EQUAL(data.ends, 1);
    BOOST_CHECK_EQUAL(G4StateManager::GetStateManager()->GetCurrentState(),
                      G4State_Idle);
  }

  BOOST_TEST_CONTEXT("end-run exception still restores Idle") {
    data = {};
    manager->beginRun();
    manager->processEvent();
    data.throwAtEnd = true;
    BOOST_CHECK_THROW(manager->endRun(), std::runtime_error);
    BOOST_CHECK(!manager->isRunActive());
    BOOST_CHECK_EQUAL(G4StateManager::GetStateManager()->GetCurrentState(),
                      G4State_Idle);
  }

  BOOST_TEST_CONTEXT("reject another thread without disturbing the run") {
    data = {};
    manager->beginRun();
    bool rejected = false;
    std::thread other([&] {
      try {
        manager->processEvent();
      } catch (const std::logic_error&) {
        rejected = true;
      }
    });
    other.join();
    BOOST_CHECK(rejected);
    BOOST_CHECK(manager->isRunActive());
    manager->processEvent();
    manager->endRun();
  }

  // The manager also closes an unfinished run when ACTS finalization is
  // skipped.
  data = {};
  manager->beginRun();
  manager->processEvent();
  manager.reset();
  BOOST_CHECK_EQUAL(data.ends, 1);
  BOOST_CHECK_EQUAL(data.processed, 1);
}

BOOST_AUTO_TEST_SUITE_END()
