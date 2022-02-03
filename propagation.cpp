#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"

#include <iostream>
using namespace ActsExamples;

using namespace Acts::UnitLiterals;

int main() {
  Sequencer::Config cfg;
  cfg.events = 500;
  cfg.numThreads = 1;

  auto field = std::make_shared<Acts::ConstantBField>(Acts::Vector3{0, 0, 2_T});

  GenericDetector detector{};
  auto [trkGeo, deco] = detector.finalize(GenericDetector::Config{}, nullptr);

  Sequencer s{cfg};

  Acts::Navigator nav{{trkGeo}};
  // Acts::EigenStepper<> stepper{field};
  Acts::AtlasStepper stepper{field};

  Acts::Propagator prop{std::move(stepper), std::move(nav)};

  auto cProp =
      std::make_shared<ActsExamples::ConcretePropagator<decltype(prop)>>(
          std::move(prop));

  auto rnd = std::make_shared<RandomNumbers>(RandomNumbers::Config{22});

  PropagationAlgorithm::Config pCfg;
  pCfg.propagatorImpl = cProp;
  pCfg.ntests = 1000;
  pCfg.sterileLogger = true;
  pCfg.covarianceTransport = true;
  pCfg.propagationStepCollection = "propagation-steps";
  pCfg.recordMaterialInteractions = false;
  pCfg.randomNumberSvc = rnd;

  s.addAlgorithm(
      std::make_shared<PropagationAlgorithm>(pCfg, Acts::Logging::INFO));

  s.run();
}
