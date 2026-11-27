// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geant4/Geant4Manager.hpp"
#include "ActsExamples/Geant4/Geant4Simulation.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
#include "ActsTests/CommonHelpers/WhiteBoardUtilities.hpp"

#include <limits>
#include <memory>
#include <stdexcept>

#include <G4PropagatorInField.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4TransportationManager.hh>
#include <G4UserEventAction.hh>

using namespace ActsExamples;

namespace {
class ObserveStepLimit final : public G4UserEventAction {
 public:
  void BeginOfEventAction(const G4Event* /*event*/) override {
    observed = G4TransportationManager::GetTransportationManager()
                   ->GetPropagatorInField()
                   ->GetLargestAcceptableStep();
    if (fail) {
      throw std::runtime_error("event failure");
    }
  }

  double observed = 0;
  bool fail = false;
};
}  // namespace

BOOST_AUTO_TEST_CASE(SequentialSessionsOnSharedHandle) {
  Geant4Simulation::Config cfg;
  TelescopeDetector::Config telescope;
  telescope.positions = {30, 60};
  telescope.stereos = {0, 0};
  cfg.detector = std::make_shared<TelescopeDetector>(telescope);
  cfg.inputParticles = "input";
  cfg.randomNumbers = std::make_shared<RandomNumbers>(RandomNumbers::Config{});
  cfg.physicsList = "MaterialPhysicsList";
  cfg.propagatorLargestAcceptableStep = 17 * Acts::UnitConstants::mm;
  cfg.magneticField =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0, 0, 1));
  std::weak_ptr<const Acts::MagneticFieldProvider> field = cfg.magneticField;
  constexpr double originalStep = 125 * CLHEP::mm;
  ObserveStepLimit* observer = nullptr;
  {
    Geant4Simulation initial(cfg);
    cfg.geant4Handle = initial.geant4Handle();
    cfg.magneticField.reset();
    observer = new ObserveStepLimit();
    cfg.geant4Handle->runManager->SetUserAction(observer);
    BOOST_REQUIRE(initial.initialize() == ProcessCode::SUCCESS);
    auto* propagator = G4TransportationManager::GetTransportationManager()
                           ->GetPropagatorInField();
    propagator->SetLargestAcceptableStep(originalStep);
    WhiteBoard first;
    ActsTests::addToWhiteBoard(cfg.inputParticles, SimParticleContainer{},
                               first);
    BOOST_REQUIRE(initial.execute({0, 0, first, 0}) == ProcessCode::SUCCESS);
    BOOST_CHECK_EQUAL(observer->observed, 17 * CLHEP::mm);
    BOOST_CHECK_EQUAL(propagator->GetLargestAcceptableStep(), originalStep);

    observer->fail = true;
    BOOST_CHECK_THROW(initial.execute({0, 1, first, 0}), std::runtime_error);
    BOOST_CHECK_EQUAL(propagator->GetLargestAcceptableStep(), originalStep);
    observer->fail = false;
    BOOST_REQUIRE(initial.initialize() == ProcessCode::SUCCESS);
    // Destruction must finalize a run even if the framework skips finalize().
  }
  BOOST_CHECK(!field.expired());

  cfg.propagatorLargestAcceptableStep = std::numeric_limits<double>::infinity();
  Geant4Simulation superseded(cfg);
  BOOST_CHECK(field.expired());
  Geant4Simulation current(cfg);
  auto invalidSimulation = cfg;
  invalidSimulation.outputParticles.clear();
  BOOST_CHECK_THROW((void)Geant4Simulation(invalidSimulation),
                    std::invalid_argument);
  Geant4MaterialRecording::Config invalidMaterial;
  invalidMaterial.detector = cfg.detector;
  invalidMaterial.randomNumbers = cfg.randomNumbers;
  invalidMaterial.inputParticles = cfg.inputParticles;
  invalidMaterial.geant4Handle = cfg.geant4Handle;
  invalidMaterial.outputMaterialTracks.clear();
  BOOST_CHECK_THROW((void)Geant4MaterialRecording(invalidMaterial),
                    std::invalid_argument);
  BOOST_CHECK_THROW(superseded.initialize(), std::logic_error);
  BOOST_REQUIRE(current.initialize() == ProcessCode::SUCCESS);
  BOOST_CHECK_THROW((void)Geant4Simulation(cfg), std::logic_error);
  WhiteBoard board;
  ActsTests::addToWhiteBoard(cfg.inputParticles, SimParticleContainer{}, board);
  BOOST_REQUIRE(current.execute({0, 0, board, 0}) == ProcessCode::SUCCESS);
  BOOST_CHECK_EQUAL(observer->observed, originalStep);
  BOOST_REQUIRE(current.finalize() == ProcessCode::SUCCESS);
}
