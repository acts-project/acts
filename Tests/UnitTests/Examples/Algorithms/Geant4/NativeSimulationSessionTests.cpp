// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/Geant4Simulation.hpp"
#include "ActsExamples/Geant4/RegionCreator.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
#include "ActsTests/CommonHelpers/WhiteBoardUtilities.hpp"

#include <chrono>
#include <functional>
#include <future>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <utility>

#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include <G4UserSteppingAction.hh>

namespace {
using namespace ActsExamples;
using namespace std::chrono_literals;

struct DiagnosticError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct Gate {
  std::promise<void> atStep;
  std::promise<void> otherAtStep;
  std::promise<void> proceed;
  std::shared_future<void> release = proceed.get_future().share();
  std::function<void()> nested;
};

class Diagnostic final : public G4UserSteppingAction {
 public:
  Diagnostic(std::shared_ptr<const Geant4::EventContext> context, Gate& gate)
      : m_context(std::move(context)), m_gate(gate) {}

  void UserSteppingAction(const G4Step* /*step*/) override {
    const auto* event =
        G4EventManager::GetEventManager()->GetConstCurrentEvent();
    const int id = event->GetEventID();
    if (m_context->current().inputParticles.begin()->particleId().particle() !=
        static_cast<unsigned int>(id)) {
      throw std::logic_error("Caller input and Geant4 event ID disagree");
    }
    if (id == 17) {
      throw DiagnosticError("Diagnostic rejected event 17");
    }
    if (id == 21) {
      m_gate.nested();
    }
    if ((id == 18 || id == 23) && !m_held) {
      m_held = true;
      (id == 18 ? m_gate.atStep : m_gate.otherAtStep).set_value();
      if (m_gate.release.wait_for(10s) != std::future_status::ready) {
        throw std::runtime_error("Timed out waiting to release event");
      }
    }
  }

 private:
  std::shared_ptr<const Geant4::EventContext> m_context;
  Gate& m_gate;
  bool m_held = false;
};

Geant4Simulation::Config makeConfig() {
  TelescopeDetector::Config telescope;
  telescope.positions = {30, 60};
  telescope.stereos = {0, 0};
  auto detector = std::make_shared<TelescopeDetector>(telescope);
  Geant4Simulation::Config cfg;
  cfg.inputParticles = "input";
  cfg.detector = detector;
  cfg.randomNumbers = std::make_shared<RandomNumbers>(RandomNumbers::Config{});
  cfg.physicsList = "MaterialPhysicsList";
  cfg.numGeant4Workers = 2;
  cfg.recordHitsOfNeutrals = true;
  Geant4::SensitiveSurfaceMapper::Config map;
  map.volumeMappings = {"Layer Logic"};
  map.candidateSurfaces = std::make_shared<Geant4::SensitiveCandidates>(
      detector->trackingGeometry(),
      Acts::getDefaultLogger("Candidates", Acts::Logging::WARNING));
  cfg.sensitiveSurfaceMapper =
      std::make_shared<Geant4::SensitiveSurfaceMapper>(map);
  return cfg;
}

std::size_t transport(Geant4Simulation& simulation, unsigned int id) {
  const auto barcode = SimBarcode{}.withVertexPrimary(1).withParticle(id);
  SimParticle particle(barcode, Acts::PdgParticle::eInvalid, 0, 0);
  particle.initialState().setDirection(0, 0, 1).setAbsoluteMomentum(1.);
  WhiteBoard board;
  ActsTests::addToWhiteBoard(simulation.config().inputParticles,
                             SimParticleContainer{particle}, board);
  if (simulation.execute({0, id, board, 0}) != ProcessCode::SUCCESS) {
    throw std::runtime_error("Geant4 execution failed");
  }
  const auto& hits = ActsTests::getFromWhiteBoard<SimHitContainer>(
      simulation.config().outputSimHits, board);
  for (const auto& hit : hits) {
    if (hit.particleId() != barcode) {
      throw std::logic_error("Hit belongs to another request");
    }
  }
  return hits.size();
}
}  // namespace

BOOST_AUTO_TEST_CASE(NativeSessionLifecycle) {
  const auto& suite = boost::unit_test::framework::master_test_suite();
  const std::string_view mode = suite.argc > 1 ? suite.argv[1] : "transport";
  Gate gate;
  auto cfg = makeConfig();
  if (mode == "empty") {
    Geant4::RegionCreator::Config region;
    region.name = "TelescopeRegion";
    region.volumes = {"Envelope Logic", "Layer Logic"};
    cfg.constructionOptions.regionCreators.push_back(
        std::make_shared<Geant4::RegionCreator>(region));
  }
  cfg.steppingActionFactories.emplace_back(
      [&](std::shared_ptr<const Geant4::EventContext> context) {
        if (mode == "startup-failure") {
          throw DiagnosticError("Diagnostic initialization failed");
        }
        return std::make_unique<Diagnostic>(std::move(context), gate);
      });
  Geant4Simulation simulation(cfg);
  gate.nested = [&] { transport(simulation, 22); };
  if (mode == "startup-failure") {
    BOOST_CHECK_THROW(simulation.initialize(), DiagnosticError);
    BOOST_CHECK_NO_THROW(simulation.finalize());
    return;
  }

  BOOST_REQUIRE(simulation.initialize() == ProcessCode::SUCCESS);
  if (mode == "empty") {
    auto* regions = G4RegionStore::GetInstance();
    const auto* region = regions->GetRegion("TelescopeRegion", false);
    BOOST_REQUIRE(region != nullptr);
    BOOST_REQUIRE_EQUAL(region->GetNumberOfRootVolumes(), 2u);
    BOOST_REQUIRE(simulation.finalize() == ProcessCode::SUCCESS);
    region = regions->GetRegion("TelescopeRegion", false);
    BOOST_CHECK(region == nullptr || region->GetNumberOfRootVolumes() == 0);
  } else {
    // Hold event 18 on one native worker. A second worker must reject event
    // 17 and recover for a subsequent event while the first is still tracking.
    auto active = std::async(std::launch::async,
                             [&] { return transport(simulation, 18); });
    BOOST_REQUIRE(gate.atStep.get_future().wait_for(10s) ==
                  std::future_status::ready);
    BOOST_CHECK_THROW(transport(simulation, 17), DiagnosticError);
    BOOST_CHECK_THROW(transport(simulation, 21), std::logic_error);
    BOOST_CHECK_EQUAL(transport(simulation, 19), 2u);

    // With both slots actively tracking, the next submission must wait for
    // admission. Closing rejects that submission while both inputs stay live.
    auto otherActive = std::async(std::launch::async,
                                  [&] { return transport(simulation, 23); });
    BOOST_REQUIRE(gate.otherAtStep.get_future().wait_for(10s) ==
                  std::future_status::ready);
    auto blocked = std::async(std::launch::async,
                              [&] { return transport(simulation, 24); });
    BOOST_CHECK(blocked.wait_for(10ms) == std::future_status::timeout);
    auto finish =
        std::async(std::launch::async, [&] { return simulation.finalize(); });
    BOOST_CHECK(blocked.wait_for(2s) == std::future_status::ready);
    BOOST_CHECK(finish.wait_for(10ms) == std::future_status::timeout);
    gate.proceed.set_value();
    BOOST_CHECK_THROW(blocked.get(), std::logic_error);
    BOOST_CHECK_EQUAL(active.get(), 2u);
    BOOST_CHECK_EQUAL(otherActive.get(), 2u);
    BOOST_REQUIRE(finish.get() == ProcessCode::SUCCESS);
  }
  BOOST_CHECK_THROW(transport(simulation, 20), std::logic_error);
  BOOST_CHECK(simulation.finalize() == ProcessCode::SUCCESS);
}
