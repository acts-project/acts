// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Digitization/HitSmearing.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

#include <boost/program_options.hpp>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int runSeedingExample(int argc, char* argv[],
                      std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addOutputOptions(desc);
  Options::addInputOptions(desc);

  // Add specific options for this geometry
  detector->addOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());

  // The geometry, material and decoration
  auto geometry = Geometry::build(vm, *detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vm));

  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Read particles (initial states) and clusters from CSV files
  auto particleReader = Options::readCsvParticleReaderConfig(vm);
  particleReader.inputStem = "particles_initial";
  particleReader.outputParticles = "particles_initial";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(particleReader, logLevel));

  // Read truth hits from CSV files
  auto simHitReaderCfg = Options::readCsvSimHitReaderConfig(vm);
  simHitReaderCfg.inputStem = "simhits";
  simHitReaderCfg.outputSimHits = "simhits";
  sequencer.addReader(
      std::make_shared<CsvSimHitReader>(simHitReaderCfg, logLevel));

  // Create smeared measurements
  HitSmearing::Config hitSmearingCfg;
  hitSmearingCfg.inputSimHits = simHitReaderCfg.outputSimHits;
  hitSmearingCfg.outputSourceLinks = "sourcelinks";
  hitSmearingCfg.outputMeasurements = "measurements";
  hitSmearingCfg.outputMeasurementParticlesMap = "measurement_particles_map";
  hitSmearingCfg.outputMeasurementSimHitsMap = "measurement_simhits_map";
  hitSmearingCfg.sigmaLoc0 = 25_um;
  hitSmearingCfg.sigmaLoc1 = 100_um;
  hitSmearingCfg.randomNumbers = rnd;
  hitSmearingCfg.trackingGeometry = tGeometry;
  sequencer.addAlgorithm(
      std::make_shared<HitSmearing>(hitSmearingCfg, logLevel));

  // Create space points
  SpacePointMaker::Config spCfg;
  spCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  spCfg.inputMeasurements = hitSmearingCfg.outputMeasurements;
  spCfg.outputSpacePoints = "spacepoints";
  spCfg.trackingGeometry = tGeometry;
  spCfg.geometrySelection = {
      // barrel pixel layers
      Acts::GeometryIdentifier().setVolume(8).setLayer(2),
      Acts::GeometryIdentifier().setVolume(8).setLayer(4),
      Acts::GeometryIdentifier().setVolume(8).setLayer(6),
      // positive endcap pixel layers
      Acts::GeometryIdentifier().setVolume(9).setLayer(2),
      Acts::GeometryIdentifier().setVolume(9).setLayer(4),
      Acts::GeometryIdentifier().setVolume(9).setLayer(6),
      Acts::GeometryIdentifier().setVolume(9).setLayer(8),
      // negative endcap pixel layers
      Acts::GeometryIdentifier().setVolume(7).setLayer(14),
      Acts::GeometryIdentifier().setVolume(7).setLayer(12),
      Acts::GeometryIdentifier().setVolume(7).setLayer(10),
      Acts::GeometryIdentifier().setVolume(7).setLayer(8),
  };
  sequencer.addAlgorithm(std::make_shared<SpacePointMaker>(spCfg, logLevel));

  // Seeding algorithm
  SeedingAlgorithm::Config seedingCfg;
  seedingCfg.inputSpacePoints = {
      spCfg.outputSpacePoints,
  };
  seedingCfg.outputSeeds = "seeds";
  seedingCfg.outputProtoTracks = "prototracks";
  seedingCfg.rMax = 200.;
  seedingCfg.deltaRMax = 60.;
  seedingCfg.collisionRegionMin = -250;
  seedingCfg.collisionRegionMax = 250.;
  seedingCfg.zMin = -2000.;
  seedingCfg.zMax = 2000.;
  seedingCfg.maxSeedsPerSpM = 1;
  seedingCfg.cotThetaMax = 7.40627;  // 2.7 eta
  seedingCfg.sigmaScattering = 50;
  seedingCfg.radLengthPerSeed = 0.1;
  seedingCfg.minPt = 500.;
  seedingCfg.bFieldInZ = 0.00199724;
  seedingCfg.beamPosX = 0;
  seedingCfg.beamPosY = 0;
  seedingCfg.impactMax = 3.;
  sequencer.addAlgorithm(
      std::make_shared<SeedingAlgorithm>(seedingCfg, logLevel));

  // Performance Writer
  TrackFinderPerformanceWriter::Config tfPerfCfg;
  tfPerfCfg.inputProtoTracks = seedingCfg.outputProtoTracks;
  tfPerfCfg.inputParticles = particleReader.outputParticles;
  tfPerfCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  tfPerfCfg.outputDir = outputDir;
  tfPerfCfg.outputFilename = "performance_seeding_trees.root";
  sequencer.addWriter(
      std::make_shared<TrackFinderPerformanceWriter>(tfPerfCfg, logLevel));

  SeedingPerformanceWriter::Config seedPerfCfg;
  seedPerfCfg.inputSeeds = seedingCfg.outputSeeds;
  seedPerfCfg.inputParticles = particleReader.outputParticles;
  seedPerfCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  seedPerfCfg.outputDir = outputDir;
  seedPerfCfg.outputFilename = "performance_seeding_hists.root";
  sequencer.addWriter(
      std::make_shared<SeedingPerformanceWriter>(seedPerfCfg, logLevel));

  return sequencer.run();
}
