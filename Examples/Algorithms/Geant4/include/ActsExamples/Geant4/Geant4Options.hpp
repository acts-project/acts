// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/Utilities/Options.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

#include <iostream>

#include "GeantinoRecording.hpp"

namespace po = boost::program_options;
using namespace Acts::UnitLiterals;

namespace FW {

namespace Options {

/// @brief ExtrapolationAlgorithm options
///
/// @tparam aopt_t Type of the options class from boost
template <typename aopt_t>
void addGeant4Options(aopt_t& opt) {
  opt.add_options()("g4-rnd-seed1",
                    po::value<unsigned int>()->default_value(287362910),
                    "The first seed of the G4 random number generation")(
      "g4-rnd-seed2", po::value<unsigned int>()->default_value(730284537),
      "The second seed of the G4 random number generation")(
      "g4-pg-nparticles", po::value<unsigned int>()->default_value(100),
      "The number of particles produced by the g4 particle gun")(
      "g4-material-tracks",
      po::value<std::string>()->default_value("geant4-material-tracks"),
      "The output collection for material tracks");
}

/// Read the Geatn4 options and @return a GeantinoRecording::Config
///
/// @tparam vmap_t is the Type of the Parameter map to be read out
///
/// @param vm is the parameter map for the options
///
/// @returns a Config object for the GeantinoRecording
template <typename vmap_t>
ActsExamples::GeantinoRecording::Config readGeantinoRecordingConfig(
    const vmap_t& vm) {
  ActsExamples::GeantinoRecording::Config gRecConfig;

  gRecConfig.tracksPerEvent =
      vm["g4-pg-nparticles"].template as<unsigned int>();
  gRecConfig.seed1 = vm["g4-rnd-seed1"].template as<unsigned int>();
  gRecConfig.seed2 = vm["g4-rnd-seed2"].template as<unsigned int>();
  gRecConfig.outputMaterialTracks =
      vm["g4-material-tracks"].template as<std::string>();

  return gRecConfig;
}

}  // namespace Options
}  // namespace FW
