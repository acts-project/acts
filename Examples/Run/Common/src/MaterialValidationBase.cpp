// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/StraightLineStepper.hpp>
#include <boost/program_options.hpp>
#include <memory>

#include "ACTFW/Detector/IBaseDetector.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Root/RootMaterialTrackWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Plugins/BField/ScalableBField.hpp"
#include "ACTFW/Propagation/PropagationAlgorithm.hpp"
#include "ACTFW/Propagation/PropagationOptions.hpp"
#include "ACTFW/Utilities/Paths.hpp"

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
template <typename bfield_t>
FW::ProcessCode setupPropagation(
    FW::Sequencer& sequencer, bfield_t bfield, po::variables_map& vm,
    std::shared_ptr<FW::RandomNumbers> randomNumberSvc,
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry) {
  // Get the log level
  auto logLevel = FW::Options::readLogLevel(vm);

  // Get a Navigator
  Acts::Navigator navigator(tGeometry);

  // Resolve the bfield map template and create the propgator
  using Stepper = Acts::EigenStepper<bfield_t>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
  Stepper stepper(std::move(bfield));
  Propagator propagator(std::move(stepper), std::move(navigator));

  // Read the propagation config and create the algorithms
  auto pAlgConfig = FW::Options::readPropagationConfig(vm, propagator);
  pAlgConfig.randomNumberSvc = randomNumberSvc;
  auto propagationAlg = std::make_shared<FW::PropagationAlgorithm<Propagator>>(
      pAlgConfig, logLevel);

  // Add the propagation algorithm
  sequencer.addAlgorithm({propagationAlg});

  return FW::ProcessCode::SUCCESS;
}

/// @brief Straight Line Propagation setup
///
/// @param sequencer The framework sequencer, Propagation algorithm to be added
/// @param vm The program options for the log file
/// @param randomNumberSvc The framework random number engine
/// @param tGeometry The TrackingGeometry object
///
/// @return a process code
FW::ProcessCode setupStraightLinePropagation(
    FW::Sequencer& sequencer, po::variables_map& vm,
    std::shared_ptr<FW::RandomNumbers> randomNumberSvc,
    std::shared_ptr<const Acts::TrackingGeometry> tGeometry) {
  // Get the log level
  auto logLevel = FW::Options::readLogLevel(vm);

  // Get a Navigator
  Acts::Navigator navigator(tGeometry);

  // Straight line stepper
  using SlStepper = Acts::StraightLineStepper;
  using Propagator = Acts::Propagator<SlStepper, Acts::Navigator>;
  // Make stepper and propagator
  SlStepper stepper;
  Propagator propagator(std::move(stepper), std::move(navigator));

  // Read the propagation config and create the algorithms
  auto pAlgConfig = FW::Options::readPropagationConfig(vm, propagator);
  pAlgConfig.randomNumberSvc = randomNumberSvc;
  auto propagationAlg = std::make_shared<FW::PropagationAlgorithm<Propagator>>(
      pAlgConfig, logLevel);

  // Add the propagation algorithm
  sequencer.addAlgorithm({propagationAlg});

  return FW::ProcessCode::SUCCESS;
}

}  // namespace

int materialValidationExample(int argc, char* argv[],
                              FW::IBaseDetector& detector) {
  // Setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addBFieldOptions(desc);
  FW::Options::addRandomNumbersOptions(desc);
  FW::Options::addPropagationOptions(desc);
  FW::Options::addOutputOptions(desc);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  FW::Sequencer sequencer(FW::Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = FW::Options::readLogLevel(vm);

  // The geometry, material and decoration
  auto geometry = FW::Geometry::build(vm, detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;

  // Create the random number engine
  auto randomNumberSvcCfg = FW::Options::readRandomNumbersConfig(vm);
  auto randomNumberSvc =
      std::make_shared<FW::RandomNumbers>(randomNumberSvcCfg);

  // Create BField service
  auto bFieldVar = FW::Options::readBField(vm);

  if (vm["prop-stepper"].template as<int>() == 0) {
    // Straight line stepper was chosen
    setupStraightLinePropagation(sequencer, vm, randomNumberSvc, tGeometry);
  } else {
    std::visit(
        [&](auto& bField) {
          using field_type =
              typename std::decay_t<decltype(bField)>::element_type;
          Acts::SharedBField<field_type> fieldMap(bField);
          setupPropagation(sequencer, fieldMap, vm, randomNumberSvc, tGeometry);
        },
        bFieldVar);
  }

  // ---------------------------------------------------------------------------------
  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();
  auto matCollection = vm["prop-material-collection"].as<std::string>();

  if (vm["output-root"].template as<bool>()) {
    // Write the propagation steps as ROOT TTree
    FW::RootMaterialTrackWriter::Config matTrackWriterRootConfig;
    matTrackWriterRootConfig.collection = matCollection;
    matTrackWriterRootConfig.filePath =
        FW::joinPaths(outputDir, matCollection + ".root");
    matTrackWriterRootConfig.storesurface = true;
    auto matTrackWriterRoot = std::make_shared<FW::RootMaterialTrackWriter>(
        matTrackWriterRootConfig, logLevel);
    sequencer.addWriter(matTrackWriterRoot);
  }

  // Initiate the run
  sequencer.run();
  // Return success code
  return 0;
}
