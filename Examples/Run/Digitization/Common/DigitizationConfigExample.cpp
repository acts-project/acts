// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Digitization/DigitizationConfigurator.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/DigitizationOptions.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <fstream>
#include <memory>

#include <boost/program_options.hpp>

#include "DigitizationInput.hpp"

using namespace ActsExamples;

int runDigitizationConfigExample(
    int argc, char* argv[],
    std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addDigitizationOptions(desc);

  auto opt = desc.add_options();
  opt("digi-compactify-output", boost::program_options::bool_switch(),
      "Try to compactify the resulting output .json file by "
      "identifying common items.");

  // Add specific options for this geometry
  detector->addOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // Get the input configuration
  auto inputConfig =
      readDigiConfigFromJson(vm["digi-config-file"].as<std::string>());

  // The geometry, material and decoration
  auto geometry = Geometry::build(vm, *detector);

  // Build a parser and visit the geometry
  ActsExamples::DigitizationConfigurator digiConfigurator;
  digiConfigurator.compactify = vm["digi-compactify-output"].as<bool>();
  digiConfigurator.inputDigiComponents = inputConfig;

  geometry.first->visitSurfaces(digiConfigurator);

  Acts::GeometryHierarchyMap<DigiComponentsConfig> outputConfig(
      digiConfigurator.outputDigiComponents);

  if (not vm["dump-digi-config"].as<std::string>().empty()) {
    writeDigiConfigToJson(outputConfig,
                          vm["dump-digi-config"].as<std::string>());
  }

  return 0;
}
