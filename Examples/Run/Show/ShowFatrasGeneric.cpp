// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file An example utility to read and print Fatras CSV files.
///
/// This examples shows how to use the framework CSV I/O and printer algorithms.

#include "ActsExamples/Detector/GenericDetectorWithOptions.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Printers/HitsPrinter.hpp"
#include "ActsExamples/Printers/ParticlesPrinter.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <memory>

#include <boost/program_options.hpp>

using namespace ActsExamples;
using boost::program_options::value;

int main(int argc, char* argv[]) {
  GenericDetectorWithOptions detector;

  // setup and parse options
  auto desc = Options::makeDefaultOptions("Read and print Fatras CSV files");
  auto opts = desc.add_options();
  Options::addSequencerOptions(desc);
  opts("input-dir", value<std::string>()->default_value(""), "");
  Options::addGeometryOptions(desc);
  detector.addOptions(desc);
  Options::addMaterialOptions(desc);
  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vars));

  // read options
  auto logLevel = Options::readLogLevel(vars);
  auto inputDir = vars["input-dir"].as<std::string>();

  // setup detector
  auto [trackingGeometry, contextDecorators] = Geometry::build(vars, detector);
  for (const auto& cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // read initial/final particles
  CsvParticleReader::Config readInitialCfg;
  readInitialCfg.inputDir = inputDir;
  readInitialCfg.inputStem = "particles_initial";
  readInitialCfg.outputParticles = "particles_initial";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(readInitialCfg, logLevel));
  CsvParticleReader::Config readFinalCfg;
  readFinalCfg.inputDir = inputDir;
  readFinalCfg.inputStem = "particles_final";
  readFinalCfg.outputParticles = "particles_final";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(readFinalCfg, logLevel));

  // read clusters/hits
  CsvPlanarClusterReader::Config readClustersCfg;
  readClustersCfg.inputDir = inputDir;
  readClustersCfg.outputClusters = "clusters";
  readClustersCfg.outputMeasurementParticlesMap = "hit_particle_map";
  readClustersCfg.outputHitIds = "hit_ids";
  readClustersCfg.outputSimHits = "simulated_hits";
  readClustersCfg.trackingGeometry = trackingGeometry;
  sequencer.addReader(std::make_shared<ActsExamples::CsvPlanarClusterReader>(
      readClustersCfg, logLevel));

  // print event data
  ParticlesPrinter::Config printInitialCfg;
  printInitialCfg.inputParticles = readInitialCfg.outputParticles;
  sequencer.addAlgorithm(
      std::make_shared<ParticlesPrinter>(printInitialCfg, logLevel));
  ParticlesPrinter::Config printFinalCfg;
  printFinalCfg.inputParticles = readFinalCfg.outputParticles;
  sequencer.addAlgorithm(
      std::make_shared<ParticlesPrinter>(printFinalCfg, logLevel));
  HitsPrinter::Config printHitsCfg;
  printHitsCfg.inputClusters = readClustersCfg.outputClusters;
  printHitsCfg.inputMeasurementParticlesMap =
      readClustersCfg.outputMeasurementParticlesMap;
  printHitsCfg.inputHitIds = readClustersCfg.outputHitIds;
  // print all hits in the container
  printHitsCfg.selectIndexLength = SIZE_MAX;
  // print all hits within a volume/layer/volume
  // printHitsCfg.selectVolume = 9;
  // printHitsCfg.selectLayer = 6;
  // printHits.selectModule = 116;
  sequencer.addAlgorithm(std::make_shared<HitsPrinter>(printHitsCfg, logLevel));

  return sequencer.run();
}
