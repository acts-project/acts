// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include <Acts/Definitions/Units.hpp>

#include <cstdlib>
#include <utility>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace ActsExamples {

namespace Options {

/// the particle gun options, the are prefixes with gp
template <typename aopt_t>
void addDD4hepOptions(aopt_t& opt) {
  opt.add_options()(
      "dd4hep-input", po::value<std::vector<std::string>>(),
      "The locations of the input DD4hep files, use 'file:foo.xml'. In case "
      "you want to read in multiple files, add the option multiple times.")(
      "dd4hep-envelopeR", po::value<double>()->default_value(1.),
      "The envelop cover in R for DD4hep volumes in mm.")(
      "dd4hep-envelopeR", po::value<double>()->default_value(1.),
      "The tolerance added to the geometrical extension in r of the "
      "layers contained to build the volume envelope around in mm.")(
      "dd4hep-envelopeZ", po::value<double>()->default_value(1.),
      "The tolerance added to the geometrical extension in z of the "
      "layers contained to build the volume envelope around in mm.")(
      "dd4hep-layerThickness", po::value<double>()->default_value(10e-10),
      "In case no surfaces (to be contained by the layer) are handed over, "
      "the layer thickness will be set to this value.")(
      "dd4hep-buildFCChh", po::value<bool>()->default_value(true),
      "If you are not building the FCChh detector please set this flag to "
      "false.")("dd4hep-loglevel", po::value<size_t>()->default_value(2),
                "The output log level of the geometry building. Please set the "
                "wished "
                "number (0 = VERBOSE, 1 = "
                "DEBUG, 2 = INFO, 3 = WARNING, 4 = ERROR, 5 = FATAL).");
}

/// read the particle gun options and return a Config file
template <typename amap_t>
ActsExamples::DD4hep::DD4hepGeometryService::Config readDD4hepConfig(
    const amap_t& vm) {
  using namespace Acts::UnitLiterals;
  ActsExamples::DD4hep::DD4hepGeometryService::Config gsConfig;
  gsConfig.logLevel =
      Acts::Logging::Level(vm["dd4hep-loglevel"].template as<size_t>());
  gsConfig.xmlFileNames =
      vm["dd4hep-input"].template as<std::vector<std::string>>();
  gsConfig.bTypePhi = Acts::equidistant;
  gsConfig.bTypeR = Acts::arbitrary;
  gsConfig.bTypeZ = Acts::equidistant;
  gsConfig.envelopeR = vm["dd4hep-envelopeR"].template as<double>() * 1_mm;
  gsConfig.envelopeZ = vm["dd4hep-envelopeZ"].template as<double>() * 1_mm;
  gsConfig.defaultLayerThickness =
      vm["dd4hep-layerThickness"].template as<double>() * 1_mm;
  if (vm["dd4hep-buildFCChh"].template as<bool>()) {
    gsConfig.sortDetectors = ActsExamples::DD4hep::sortFCChhDetElements;
  }
  return gsConfig;
}
}  // namespace Options
}  // namespace ActsExamples
