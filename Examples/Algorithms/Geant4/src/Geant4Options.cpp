// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/Geant4Options.hpp"

#include <boost/program_options.hpp>

#include <string>

void ActsExamples::Options::addGeant4Options(
    ActsExamples::Options::Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("g4-rnd-seed1", value<unsigned int>()->default_value(287362910),
      "The first seed of the G4 random number generation");
  opt("g4-rnd-seed2", value<unsigned int>()->default_value(730284537),
      "The second seed of the G4 random number generation");
  opt("g4-pg-nparticles", value<unsigned int>()->default_value(100),
      "The number of particles produced by the g4 particle gun");
  opt("g4-material-tracks",
      value<std::string>()->default_value("geant4-material-tracks"),
      "The output collection for material tracks");
}

ActsExamples::GeantinoRecording::Config
ActsExamples::Options::readGeantinoRecordingConfig(
    const ActsExamples::Options::Variables& variables) {
  ActsExamples::GeantinoRecording::Config gRecConfig;

  gRecConfig.tracksPerEvent = variables["g4-pg-nparticles"].as<unsigned int>();
  gRecConfig.seed1 = variables["g4-rnd-seed1"].as<unsigned int>();
  gRecConfig.seed2 = variables["g4-rnd-seed2"].as<unsigned int>();
  gRecConfig.outputMaterialTracks =
      variables["g4-material-tracks"].as<std::string>();

  return gRecConfig;
}
