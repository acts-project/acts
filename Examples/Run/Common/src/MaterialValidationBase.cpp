// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/DefaultExtension.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsExamples/Propagation/PropagationOptions.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/StraightLineStepper.hpp>

#include <memory>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace {
/// @brief Propagation setup
///
/// @tparam bfield_t Type of the magnetic field
///
/// @param sequencer The framework sequencer, Propagation algorithm to be added
/// @param bfield The bfield object needed for the Stepper & propagagor
/// @param vm The program options for the log file
/// @param randomNumberSvc The framework random number engine
/// @param tGeometry The TrackingGeometry object
///
/// @return a process code
ActsExamples::ProcessCode setupPropagation(
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const Acts::MagneticFieldProvider> bfield,
    po::variables_map& vm,
    std::shared_ptr<ActsExamples::RandomNumbers> randomNumberSvc,
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry) {
  // Get the log level
  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // Get a Navigator
  Acts::Navigator::Config cfg;
  cfg.trackingGeometry = tGeometry;
  cfg.resolvePassive = true;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);

  // Resolve the bfield map template and create the propgator
  using Stepper = Acts::EigenStepper<
      Acts::StepperExtensionList<Acts::DefaultExtension,
                                 Acts::DenseEnvironmentExtension>,
      Acts::detail::HighestValidAuctioneer>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
  Stepper stepper(std::move(bfield));
  Propagator propagator(std::move(stepper), std::move(navigator));

  // Read the propagation config and create the algorithms
  auto pAlgConfig = ActsExamples::Options::readPropagationConfig(vm);
  pAlgConfig.randomNumberSvc = randomNumberSvc;
  pAlgConfig.recordMaterialInteractions = true;

  pAlgConfig.propagatorImpl =
      std::make_shared<ActsExamples::ConcretePropagator<Propagator>>(
          std::move(propagator));

  auto propagationAlg = std::make_shared<ActsExamples::PropagationAlgorithm>(
      pAlgConfig, logLevel);

  // Add the propagation algorithm
  sequencer.addAlgorithm({propagationAlg});

  return ActsExamples::ProcessCode::SUCCESS;
}

/// @brief Straight Line Propagation setup
///
/// @param sequencer The framework sequencer, Propagation algorithm to be added
/// @param vm The program options for the log file
/// @param randomNumberSvc The framework random number engine
/// @param tGeometry The TrackingGeometry object
///
/// @return a process code
ActsExamples::ProcessCode setupStraightLinePropagation(
    ActsExamples::Sequencer& sequencer, po::variables_map& vm,
    std::shared_ptr<ActsExamples::RandomNumbers> randomNumberSvc,
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry) {
  // Get the log level
  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // Get a Navigator
  Acts::Navigator navigator({tGeometry});

  // Straight line stepper
  using SlStepper = Acts::StraightLineStepper;
  using Propagator = Acts::Propagator<SlStepper, Acts::Navigator>;
  // Make stepper and propagator
  SlStepper stepper;
  Propagator propagator(stepper, std::move(navigator));

  // Read the propagation config and create the algorithms
  auto pAlgConfig = ActsExamples::Options::readPropagationConfig(vm);

  pAlgConfig.randomNumberSvc = randomNumberSvc;
  pAlgConfig.propagatorImpl =
      std::make_shared<ActsExamples::ConcretePropagator<Propagator>>(
          std::move(propagator));
  auto propagationAlg = std::make_shared<ActsExamples::PropagationAlgorithm>(
      pAlgConfig, logLevel);

  // Add the propagation algorithm
  sequencer.addAlgorithm({propagationAlg});

  return ActsExamples::ProcessCode::SUCCESS;
}

}  // namespace

int materialValidationExample(int argc, char* argv[],
                              ActsExamples::IBaseDetector& detector) {
  // Setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addGeometryOptions(desc);
  ActsExamples::Options::addMaterialOptions(desc);
  ActsExamples::Options::addMagneticFieldOptions(desc);
  ActsExamples::Options::addRandomNumbersOptions(desc);
  ActsExamples::Options::addPropagationOptions(desc);
  ActsExamples::Options::addOutputOptions(desc,
                                          ActsExamples::OutputFormat::Root);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

  // Now read the standard options
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

  // Create BField service
  ActsExamples::Options::setupMagneticFieldServices(vm, sequencer);
  auto bField = ActsExamples::Options::readMagneticField(vm);

  if (vm["prop-stepper"].template as<int>() == 0) {
    // Straight line stepper was chosen
    setupStraightLinePropagation(sequencer, vm, randomNumberSvc, tGeometry);
  } else {
    setupPropagation(sequencer, bField, vm, randomNumberSvc, tGeometry);
  }

  // ---------------------------------------------------------------------------------
  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();
  auto matCollection = vm["prop-material-collection"].as<std::string>();

  if (vm["output-root"].template as<bool>()) {
    // Write the propagation steps as ROOT TTree
    ActsExamples::RootMaterialTrackWriter::Config matTrackWriterRootConfig;
    matTrackWriterRootConfig.collection = matCollection;
    matTrackWriterRootConfig.filePath =
        ActsExamples::joinPaths(outputDir, matCollection + ".root");
    matTrackWriterRootConfig.storeSurface = true;
    matTrackWriterRootConfig.storeVolume = true;
    auto matTrackWriterRoot =
        std::make_shared<ActsExamples::RootMaterialTrackWriter>(
            matTrackWriterRootConfig, logLevel);
    sequencer.addWriter(matTrackWriterRoot);
  }

  // Initiate the run
  sequencer.run();
  // Return success code
  return 0;
}
