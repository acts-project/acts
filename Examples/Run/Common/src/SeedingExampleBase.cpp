// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsExamples/Seeding/SeedingAlgorithm.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Geometry/TrackingGeometry.hpp>

#include <memory>

#include <boost/program_options.hpp>

int seedingExample(int argc, char* argv[],
                   ActsExamples::IBaseDetector& detector) {
  // Setup and parse options

  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addGeometryOptions(desc);
  ActsExamples::Options::addMaterialOptions(desc);
  ActsExamples::Options::addOutputOptions(desc);
  ActsExamples::Options::addInputOptions(desc);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // The geometry, material and decoration
  auto geometry = ActsExamples::Geometry::build(vm, detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Read particles (initial states) and clusters from CSV files
  auto particleReader = ActsExamples::Options::readCsvParticleReaderConfig(vm);
  particleReader.inputStem = "particles_initial";
  particleReader.outputParticles = "particles_initial";
  sequencer.addReader(std::make_shared<ActsExamples::CsvParticleReader>(
      particleReader, logLevel));

  // Read clusters from CSV files
  auto clusterReaderCfg =
      ActsExamples::Options::readCsvPlanarClusterReaderConfig(vm);
  clusterReaderCfg.trackingGeometry = tGeometry;
  clusterReaderCfg.outputClusters = "clusters";
  clusterReaderCfg.outputHitIds = "hit_ids";
  clusterReaderCfg.outputHitParticlesMap = "hit_particles_map";
  clusterReaderCfg.outputSimulatedHits = "hits";
  sequencer.addReader(std::make_shared<ActsExamples::CsvPlanarClusterReader>(
      clusterReaderCfg, logLevel));

  // Seeding algorithm
  ActsExamples::SeedingAlgorithm::Config seeding;
  seeding.outputSeeds = "seeds";
  seeding.outputProtoTracks = "protoTracks";
  seeding.inputHitParticlesMap = clusterReaderCfg.outputHitParticlesMap;
  seeding.inputClusters = clusterReaderCfg.outputClusters;
  seeding.inputParticles = particleReader.outputParticles;
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::SeedingAlgorithm>(seeding, logLevel));

  // Performance Writer
  ActsExamples::SeedingPerformanceWriter::Config seedPerfCfg;
  // seedPerfCfg.inputSeeds = seeding.outputSeeds;
  seedPerfCfg.inputSeeds = "seeds";
  seedPerfCfg.inputProtoTracks = seeding.outputProtoTracks;
  seedPerfCfg.inputParticles = particleReader.outputParticles;
  seedPerfCfg.inputClusters = clusterReaderCfg.outputClusters;
  seedPerfCfg.inputHitParticlesMap = clusterReaderCfg.outputHitParticlesMap;
  seedPerfCfg.outputFilename = "performance.root";
  sequencer.addWriter(std::make_shared<ActsExamples::SeedingPerformanceWriter>(
      seedPerfCfg, logLevel));

  return sequencer.run();
}
