// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/AtlasStepper.hpp>
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
#include "ACTFW/Io/Root/RootPropagationStepsWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Plugins/BField/ScalableBField.hpp"
#include "ACTFW/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ACTFW/Propagation/PropagationAlgorithm.hpp"
#include "ACTFW/Propagation/PropagationOptions.hpp"
#include "ACTFW/Utilities/Paths.hpp"

int propagationExample(int argc, char* argv[], FW::IBaseDetector& detector) {
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
  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Create the random number engine
  auto randomNumberSvcCfg = FW::Options::readRandomNumbersConfig(vm);
  auto randomNumberSvc =
      std::make_shared<FW::RandomNumbers>(randomNumberSvcCfg);

  // Create BField service
  auto bFieldVar = FW::Options::readBField(vm);
  // auto field2D = std::get<std::shared_ptr<InterpolatedBFieldMap2D>>(bField);
  // auto field3D = std::get<std::shared_ptr<InterpolatedBFieldMap3D>>(bField);

  // Get a Navigator
  Acts::Navigator navigator(tGeometry);

  std::visit(
      [&](auto& bField) {
        // Resolve the bfield map and create the propgator
        using field_type =
            typename std::decay_t<decltype(bField)>::element_type;
        Acts::SharedBField<field_type> fieldMap(bField);

        using field_map_type = decltype(fieldMap);

        std::optional<std::variant<Acts::EigenStepper<field_map_type>,
                                   Acts::AtlasStepper<field_map_type>,
                                   Acts::StraightLineStepper>>
            var_stepper;

        // translate option to variant
        if (vm["prop-stepper"].template as<int>() == 0) {
          var_stepper = Acts::StraightLineStepper{};
        } else if (vm["prop-stepper"].template as<int>() == 1) {
          var_stepper = Acts::EigenStepper<field_map_type>{std::move(fieldMap)};
        } else if (vm["prop-stepper"].template as<int>() == 2) {
          var_stepper = Acts::AtlasStepper<field_map_type>{std::move(fieldMap)};
        }

        // resolve stepper, setup propagator
        std::visit(
            [&](auto& stepper) {
              using Stepper = std::decay_t<decltype(stepper)>;
              using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
              Propagator propagator(std::move(stepper), std::move(navigator));

              // Read the propagation config and create the algorithms
              auto pAlgConfig =
                  FW::Options::readPropagationConfig(vm, propagator);
              pAlgConfig.randomNumberSvc = randomNumberSvc;
              sequencer.addAlgorithm(
                  std::make_shared<FW::PropagationAlgorithm<Propagator>>(
                      pAlgConfig, logLevel));
            },
            *var_stepper);
      },
      bFieldVar);

  // ---------------------------------------------------------------------------------
  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();
  auto psCollection = vm["prop-step-collection"].as<std::string>();

  if (vm["output-root"].template as<bool>()) {
    // Write the propagation steps as ROOT TTree
    FW::RootPropagationStepsWriter::Config pstepWriterRootConfig;
    pstepWriterRootConfig.collection = psCollection;
    pstepWriterRootConfig.filePath =
        FW::joinPaths(outputDir, psCollection + ".root");
    sequencer.addWriter(std::make_shared<FW::RootPropagationStepsWriter>(
        pstepWriterRootConfig));
  }

  if (vm["output-obj"].template as<bool>()) {
    using PropagationSteps = Acts::detail::Step;
    using ObjPropagationStepsWriter =
        FW::Obj::ObjPropagationStepsWriter<PropagationSteps>;

    // Write the propagation steps as Obj TTree
    ObjPropagationStepsWriter::Config pstepWriterObjConfig;
    pstepWriterObjConfig.collection = psCollection;
    pstepWriterObjConfig.outputDir = outputDir;
    sequencer.addWriter(
        std::make_shared<ObjPropagationStepsWriter>(pstepWriterObjConfig));
  }

  return sequencer.run();
}
