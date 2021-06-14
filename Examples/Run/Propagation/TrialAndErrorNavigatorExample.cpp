// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/GuidedNavigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/TrialAndErrorSurfaceProvider.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvPropagationStepsWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsExamples/Propagation/PropagationOptions.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

int main(int argc, char** argv) {
  GenericDetector detector;

  auto main_logger = Acts::getDefaultLogger("Main", Acts::Logging::INFO);
  ACTS_LOCAL_LOGGER(std::move(main_logger));

  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addGeometryOptions(desc);
  ActsExamples::Options::addMaterialOptions(desc);
  ActsExamples::Options::addMagneticFieldOptions(desc);
  ActsExamples::Options::addRandomNumbersOptions(desc);
  ActsExamples::Options::addPropagationOptions(desc);
  ActsExamples::Options::addOutputOptions(
      desc, ActsExamples::OutputFormat::Csv | ActsExamples::OutputFormat::Obj);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // Sequencer
  const auto sequencer_config = ActsExamples::Options::readSequencerConfig(vm);
  ActsExamples::Sequencer sequencer(sequencer_config);

  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // The geometry, material and decoration
  auto geometry = ActsExamples::Geometry::build(vm, detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Create the random number engine
  auto randomNumberSvcCfg = ActsExamples::Options::readRandomNumbersConfig(vm);
  auto randomNumberSvc =
      std::make_shared<ActsExamples::RandomNumbers>(randomNumberSvcCfg);

  // Magnetic Field
  auto bFieldVar = ActsExamples::Options::readMagneticField(vm);

  // Create BField service
  ActsExamples::Options::setupMagneticFieldServices(vm, sequencer);
  auto bField = ActsExamples::Options::readMagneticField(vm);

  auto setupPropagator = [&](auto&& stepper, auto&& navigator) {
    using Stepper = std::decay_t<decltype(stepper)>;
    using Propagator =
        Acts::Propagator<Stepper, std::decay_t<decltype(navigator)>>;

    Propagator propagator(std::move(stepper), std::move(navigator));

    // Read the propagation config and create the algorithms
    auto pAlgConfig =
        ActsExamples::Options::readPropagationConfig(vm, propagator);
    pAlgConfig.randomNumberSvc = randomNumberSvc;
    sequencer.addAlgorithm(
        std::make_shared<ActsExamples::PropagationAlgorithm<Propagator>>(
            pAlgConfig, logLevel));
  };

  Acts::GuidedNavigator<Acts::TrialAndErrorSurfaceProvider> navigator(
      Acts::TrialAndErrorSurfaceProvider(*tGeometry), tGeometry);

  setupPropagator(Acts::EigenStepper<>{std::move(bField)},
                  std::move(navigator));

  // Output
  std::string outputDir = vm["output-dir"].as<std::string>();
  auto psCollection = vm["prop-step-collection"].as<std::string>();

  // Csv Writer
  if (vm["output-csv"].template as<bool>()) {
    using Writer = ActsExamples::CsvPropagationStepsWriter;

    Writer::Config config;
    config.collection = psCollection;
    config.outputDir = outputDir;

    sequencer.addWriter(std::make_shared<Writer>(config));
  }

  // Obj Writer
  if (vm["output-obj"].template as<bool>()) {
    using PropagationSteps = Acts::detail::Step;
    using ObjPropagationStepsWriter =
        ActsExamples::ObjPropagationStepsWriter<PropagationSteps>;

    // Write the propagation steps as Obj TTree
    ObjPropagationStepsWriter::Config pstepWriterObjConfig;
    pstepWriterObjConfig.collection = psCollection;
    pstepWriterObjConfig.outputDir = outputDir;
    sequencer.addWriter(
        std::make_shared<ObjPropagationStepsWriter>(pstepWriterObjConfig));
  }

  // Run sequencer
  return sequencer.run();
}
