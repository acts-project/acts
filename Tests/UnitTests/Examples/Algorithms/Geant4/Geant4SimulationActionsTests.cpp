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
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
#include "ActsTests/CommonHelpers/WhiteBoardUtilities.hpp"

#include <memory>
#include <utility>
#include <vector>

#include <G4Step.hh>
#include <G4UserSteppingAction.hh>

using namespace ActsExamples;

namespace {

class DiagnosticAction final : public G4UserSteppingAction {
 public:
  DiagnosticAction(std::shared_ptr<const Geant4::EventContext> context,
                   int index, std::vector<int>& calls)
      : m_context(std::move(context)), m_index(index), m_calls(calls) {}

  void UserSteppingAction(const G4Step* step) override {
    m_calls.push_back(m_index);
    BOOST_REQUIRE(fpSteppingManager != nullptr);
    const auto& event = m_context->current();
    const auto& track = *step->GetTrack();
    BOOST_CHECK(event.trackIdMapping.contains(track.GetTrackID()));
    // Each geantino crosses a whole sensor in a single step. Its hit must
    // already be recorded when the diagnostic receives that same step.
    if (track.GetVolume()->GetName().starts_with(
            Geant4::SensitiveSurfaceMapper::mappingPrefix)) {
      ++m_sensorSteps;
      BOOST_CHECK_EQUAL(event.result.hits.size(), m_sensorSteps);
    }
  }

 private:
  std::shared_ptr<const Geant4::EventContext> m_context;
  int m_index;
  std::vector<int>& m_calls;
  std::size_t m_sensorSteps = 0;
};

}  // namespace

BOOST_AUTO_TEST_CASE(OrderedDiagnosticsAndIndependentEventOutputs) {
  TelescopeDetector::Config telescope;
  telescope.positions = {30, 60};
  telescope.stereos = {0, 0};
  auto detector = std::make_shared<TelescopeDetector>(telescope);

  Geant4Simulation::Config cfg;
  cfg.inputParticles = "input";
  cfg.detector = detector;
  cfg.randomNumbers = std::make_shared<RandomNumbers>(RandomNumbers::Config{});
  cfg.physicsList = "MaterialPhysicsList";
  cfg.recordHitsOfNeutrals = true;
  Geant4::SensitiveSurfaceMapper::Config mapping;
  mapping.volumeMappings = {"Layer Logic"};
  mapping.candidateSurfaces = std::make_shared<Geant4::SensitiveCandidates>(
      detector->trackingGeometry(),
      Acts::getDefaultLogger("Candidates", Acts::Logging::WARNING));
  cfg.sensitiveSurfaceMapper =
      std::make_shared<Geant4::SensitiveSurfaceMapper>(mapping);

  // Invalid configuration must be rejected before acquiring the single
  // default Geant4 handle. The valid configuration below then acquires it.
  auto invalid = cfg;
  invalid.steppingActionFactories.emplace_back();
  BOOST_CHECK_THROW(Geant4Simulation{invalid}, std::invalid_argument);

  std::vector<int> constructed;
  std::vector<int> calls;
  int attempts = 0;
  for (int i : {1, 2}) {
    cfg.steppingActionFactories.emplace_back(
        [&, i](std::shared_ptr<const Geant4::EventContext> context)
            -> std::unique_ptr<G4UserSteppingAction> {
          constructed.push_back(i);
          // The first diagnostic has been constructed when the second fails.
          if (i == 2) {
            ++attempts;
            if (attempts == 1) {
              throw std::runtime_error("diagnostic unavailable");
            }
            if (attempts == 2) {
              return nullptr;
            }
          }
          return std::make_unique<DiagnosticAction>(std::move(context), i,
                                                    calls);
        });
  }
  Geant4Simulation simulation(cfg);
  BOOST_CHECK(constructed.empty());
  BOOST_CHECK_THROW(simulation.initialize(), std::runtime_error);
  BOOST_CHECK_THROW(simulation.initialize(), std::invalid_argument);
  BOOST_REQUIRE(simulation.initialize() == ProcessCode::SUCCESS);
  const std::vector<int> expectedOrder{1, 2, 1, 2, 1, 2};
  BOOST_CHECK(constructed == expectedOrder);

  const auto barcode = SimBarcode().withVertexPrimary(1).withParticle(1);
  SimParticle particle(barcode, Acts::PdgParticle::eInvalid, 0, 0);
  particle.initialState().setDirection(0, 0, 1).setAbsoluteMomentum(1);
  WhiteBoard first;
  ActsTests::addToWhiteBoard(cfg.inputParticles, SimParticleContainer{particle},
                             first);
  BOOST_REQUIRE(simulation.execute({0, 0, first, 0}) == ProcessCode::SUCCESS);
  BOOST_REQUIRE(!calls.empty());
  BOOST_REQUIRE_EQUAL(calls.size() % 2, 0u);
  for (std::size_t i = 0; i < calls.size(); i += 2) {
    BOOST_CHECK_EQUAL(calls[i], 1);
    BOOST_CHECK_EQUAL(calls[i + 1], 2);
  }
  const auto firstEventCalls = calls.size();

  WhiteBoard second;
  ActsTests::addToWhiteBoard(cfg.inputParticles, SimParticleContainer{},
                             second);
  BOOST_REQUIRE(simulation.execute({0, 1, second, 0}) == ProcessCode::SUCCESS);
  BOOST_CHECK_EQUAL(calls.size(), firstEventCalls);
  BOOST_CHECK(constructed == expectedOrder);
  BOOST_REQUIRE(simulation.finalize() == ProcessCode::SUCCESS);

  const auto particles = ActsTests::getFromWhiteBoard<SimParticleContainer>(
      cfg.outputParticles, first);
  BOOST_REQUIRE_EQUAL(particles.size(), 1u);
  BOOST_CHECK(particles.begin()->particleId() == barcode);
  const auto hits =
      ActsTests::getFromWhiteBoard<SimHitContainer>(cfg.outputSimHits, first);
  BOOST_REQUIRE_EQUAL(hits.size(), telescope.positions.size());
  for (const auto& hit : hits) {
    BOOST_CHECK(hit.particleId() == barcode);
  }
  BOOST_CHECK(ActsTests::getFromWhiteBoard<SimParticleContainer>(
                  cfg.outputParticles, second)
                  .empty());
  BOOST_CHECK(
      ActsTests::getFromWhiteBoard<SimHitContainer>(cfg.outputSimHits, second)
          .empty());
}
