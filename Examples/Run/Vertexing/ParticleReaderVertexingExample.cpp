// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file Find vertices using truth particle information as input
///
/// Reads truth particles from TrackMl files and use the truth information
/// to generate smeared track parameters. Use this pseudo-reconstructed
/// tracks as the input to the vertex finder.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"
#include "ActsExamples/Reconstruction/ReconstructionBase.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/Vertexing/IterativeVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/VertexingOptions.hpp"

#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  ParticleSelector::addOptions(desc);
  Options::addVertexingOptions(desc);
  Options::addInputOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  Options::addParticleSmearingOptions(desc);
  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // basic setup
  auto logLevel = Options::readLogLevel(vars);
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));
  Sequencer sequencer(Options::readSequencerConfig(vars));

  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vars);

  // setup particle reader generator
  CsvParticleReader::Config readParticles =
      Options::readCsvParticleReaderConfig(vars);
  readParticles.inputStem = "particles";
  readParticles.outputParticles = "particles";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(readParticles, logLevel));

  // pre-select particles
  ParticleSelector::Config selectParticles = ParticleSelector::readConfig(vars);
  selectParticles.inputParticles = readParticles.outputParticles;
  selectParticles.outputParticles = "particles_selected";
  // smearing only works with charge particles for now
  selectParticles.removeNeutral = true;
  selectParticles.absEtaMax = vars["vertexing-eta-max"].as<double>();
  selectParticles.rhoMax = vars["vertexing-rho-max"].as<double>() * 1_mm;
  selectParticles.ptMin = vars["vertexing-pt-min"].as<double>() * 1_MeV;
  sequencer.addAlgorithm(
      std::make_shared<ParticleSelector>(selectParticles, logLevel));

  // Run the particle smearing
  auto particleSmearingCfg = setupParticleSmearing(
      vars, sequencer, rnd, selectParticles.outputParticles);

  // print input track parameters
  TrackParametersPrinter::Config printTracks;
  printTracks.inputTrackParameters = particleSmearingCfg.outputTrackParameters;
  sequencer.addAlgorithm(
      std::make_shared<TrackParametersPrinter>(printTracks, logLevel));

  // find vertices
  IterativeVertexFinderAlgorithm::Config findVertices;
  findVertices.bField = magneticField;
  findVertices.inputTrackParameters = particleSmearingCfg.outputTrackParameters;
  findVertices.outputProtoVertices = "protovertices";
  sequencer.addAlgorithm(
      std::make_shared<IterativeVertexFinderAlgorithm>(findVertices, logLevel));

  return sequencer.run();
}
