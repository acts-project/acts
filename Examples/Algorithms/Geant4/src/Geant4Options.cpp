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
  opt("g4-vertex-posX", value<double>()->default_value(0.0),
      "The X position of the geantino vertex");
  opt("g4-vertex-posY", value<double>()->default_value(0.0),
      "The Y position of the geantino vertex");
  opt("g4-vertex-posZ", value<double>()->default_value(0.0),
      "The Z position of the geantino vertex");
  opt("g4-vertex-sigmaX", value<double>()->default_value(0.0),
      "The X spread of the geantino vertex");
  opt("g4-vertex-sigmaY", value<double>()->default_value(0.0),
      "The Y spread of the geantino vertex");
  opt("g4-vertex-sigmaZ", value<double>()->default_value(0.0),
      "The Z spread of the geantino vertex");
  opt("g4-phi-range",
      value<read_range>()->multitoken()->default_value({-M_PI, M_PI}),
      "Azimutal angle phi range for the geantino");
  opt("g4-eta-range",
      value<read_range>()->multitoken()->default_value({-5., 5.}),
      "Pseudorapidity eta range for the geantino");
  opt("g4-sampling-variable", value<std::string>()->default_value("theta"),
      "Variable from which the particle generation is uniform. Can be theta "
      "or eta");
}

ActsExamples::GeantinoRecording::Config
ActsExamples::Options::readGeantinoRecordingConfig(
    const ActsExamples::Options::Variables& variables) {
  ActsExamples::GeantinoRecording::Config gRecConfig;

  gRecConfig.tracksPerEvent = variables["g4-pg-nparticles"].as<unsigned int>();
  gRecConfig.generationConfig.randomSeed1 =
      variables["g4-rnd-seed1"].as<unsigned int>();
  gRecConfig.generationConfig.randomSeed2 =
      variables["g4-rnd-seed2"].as<unsigned int>();
  gRecConfig.outputMaterialTracks =
      variables["g4-material-tracks"].as<std::string>();

  gRecConfig.generationConfig.vertexPosX =
      variables["g4-vertex-posX"].as<double>();
  gRecConfig.generationConfig.vertexPosY =
      variables["g4-vertex-posY"].as<double>();
  gRecConfig.generationConfig.vertexPosZ =
      variables["g4-vertex-posZ"].as<double>();
  gRecConfig.generationConfig.vertexSigmaX =
      variables["g4-vertex-sigmaX"].as<double>();
  gRecConfig.generationConfig.vertexSigmaY =
      variables["g4-vertex-sigmaY"].as<double>();
  gRecConfig.generationConfig.vertexSigmaZ =
      variables["g4-vertex-sigmaZ"].as<double>();

  read_range iphir = variables["g4-phi-range"].template as<read_range>();
  read_range ietar = variables["g4-eta-range"].template as<read_range>();

  gRecConfig.generationConfig.phiRange = {iphir[0], iphir[1]};
  gRecConfig.generationConfig.etaRange = {ietar[0], ietar[1]};

  gRecConfig.generationConfig.samplingVariable =
      variables["g4-sampling-variable"].as<std::string>();

  return gRecConfig;
}
