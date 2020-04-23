// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <memory>

#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/GenericDetector/GenericDetector.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Csv/CsvOptionsReader.hpp"
#include "ACTFW/Io/Csv/CsvParticleReader.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Printers/PrintHits.hpp"
#include "ACTFW/Utilities/Options.hpp"

int main(int argc, char* argv[]) {
  GenericDetector detector;

  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addInputOptions(desc);
  detector.addOptions(desc);

  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  FW::Sequencer sequencer(FW::Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = FW::Options::readLogLevel(vm);
  auto inputDir = vm["input-dir"].as<std::string>();

  // Setup detector geometry
  auto geometry = FW::Geometry::build(vm, detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }

  // Read particles from CSV files
  auto particleReaderCfg = FW::Options::readCsvParticleReaderConfig(vm);
  particleReaderCfg.outputParticles = "particles";
  sequencer.addReader(
      std::make_shared<FW::CsvParticleReader>(particleReaderCfg, logLevel));

  // Read clusters from CSV files
  auto clusterReaderCfg = FW::Options::readCsvPlanarClusterReaderConfig(vm);
  clusterReaderCfg.trackingGeometry = trackingGeometry;
  clusterReaderCfg.outputClusters = "clusters";
  clusterReaderCfg.outputHitParticlesMap = "hit_particle_map";
  clusterReaderCfg.outputHitIds = "hit_ids";
  sequencer.addReader(
      std::make_shared<FW::CsvPlanarClusterReader>(clusterReaderCfg, logLevel));

  // Print some information as crosscheck
  FW::PrintHits::Config printCfg;
  printCfg.inputClusters = clusterReaderCfg.outputClusters;
  printCfg.inputHitParticlesMap = clusterReaderCfg.outputHitParticlesMap;
  printCfg.inputHitIds = clusterReaderCfg.outputHitIds;
  // the following print selections work for the original author.
  // you probably need to adapt them to your data.
  printCfg.hitIdStart = 10224;
  printCfg.hitIdLength = 8;
  printCfg.volumeId = 13;
  printCfg.layerId = 4;
  printCfg.moduleId = 116;
  sequencer.addAlgorithm(std::make_shared<FW::PrintHits>(printCfg, logLevel));

  return sequencer.run();
}
