// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Io/Csv/CsvSpacePointReader.hpp"

#include <memory>

#include <boost/program_options.hpp>

#include "RecInput.hpp"

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  
  std::cout << "I am running something..." << std::endl;
  
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  Options::addInputOptions(desc);
  Options::addMagneticFieldOptions(desc);
  
  // Add specific options for this geometry
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  
  for (const auto& it : vm) {
    std::cout << it.first.c_str() << " ";
    auto& value = it.second.value();
    if (auto v = boost::any_cast<std::string>(&value))
      std::cout << *v;
    else if (auto a = boost::any_cast<float>(&value))
      std::cout << *a;
    else if (auto b = boost::any_cast<double>(&value))
      std::cout << *b;
    else if (auto c = boost::any_cast<int>(&value))
      std::cout << *c;
    else if (auto d = boost::any_cast<bool>(&value))
      std::cout << *d;
    else if (auto e = boost::any_cast<size_t>(&value))
      std::cout << *e;
    else 
      std::cout << "wrong cast...";
    std::cout << std::endl;
  }
  
  Sequencer sequencer(Options::readSequencerConfig(vm));  
  
  // Now read the standard options
  //auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
 
  // Setup the magnetic field
  Options::setupMagneticFieldServices(vm, sequencer);
  auto magneticField = Options::readMagneticField(vm);

  //    Read the space points and build the container
  auto spReaderCfg = setupSpacePointReading(vm, sequencer);
  //spReaderCfg.outputSpacePoints = {"PixelSpacePoints", "StripSpacePoints", "OverlapSpacePoints"};
/*
  // Seeding algorithm
  SeedingAlgorithm::Config seedingCfg;
  seedingCfg.inputSpacePoints = {"PixelSpacePoints"};
  seedingCfg.outputSeeds = "PixelSeeds";
  seedingCfg.outputProtoTracks = "prototracks";
  seedingCfg.rMax = 320.;
  seedingCfg.deltaRMax = (320.-40.);
  seedingCfg.collisionRegionMin = -200;
  seedingCfg.collisionRegionMax = 200.;
  seedingCfg.zMin = -2700.;
  seedingCfg.zMax = 2700.;
  seedingCfg.maxSeedsPerSpM = 5;
  seedingCfg.cotThetaMax = 1.32582;  // 4.0 eta
  seedingCfg.sigmaScattering = 50;
  seedingCfg.radLengthPerSeed = 0.1;
  seedingCfg.minPt = 900.;
  seedingCfg.bFieldInZ = 0.00199724;
  seedingCfg.beamPosX = 0;
  seedingCfg.beamPosY = 0;
  seedingCfg.impactMax = 2.;
  seedingCfg.numberOfPhiBins = 3;
  seedingCfg.zBinEdges = {-2700, -2500., -1400., -925., -450., -250., 
                            250., 450, 925., 1400., 2500, 2700};
  
  sequencer.addAlgorithm(
      std::make_shared<SeedingAlgorithm>(seedingCfg, logLevel));

  return sequencer.run();
*/  
//   // Algorithm estimating track parameter from seed
//   TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
//   paramsEstimationCfg.inputSeeds = seedingCfg.outputSeeds;
//   paramsEstimationCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
//   paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
//   paramsEstimationCfg.outputTrackParametersSeedMap = "estimatedparams_seed_map";
//   paramsEstimationCfg.trackingGeometry = tGeometry;
//   paramsEstimationCfg.magneticField = magneticField;
//   sequencer.addAlgorithm(std::make_shared<TrackParamsEstimationAlgorithm>(
//       paramsEstimationCfg, logLevel));
// 
//   // Seeding performance Writers
//   TrackFinderPerformanceWriter::Config tfPerfCfg;
//   tfPerfCfg.inputProtoTracks = seedingCfg.outputProtoTracks;
//   tfPerfCfg.inputParticles = inputParticles;
//   tfPerfCfg.inputMeasurementParticlesMap =
//       hitSmearingCfg.outputMeasurementParticlesMap;
//   tfPerfCfg.outputDir = outputDir;
//   tfPerfCfg.outputFilename = "performance_seeding_trees.root";
//   sequencer.addWriter(
//       std::make_shared<TrackFinderPerformanceWriter>(tfPerfCfg, logLevel));
// 
//   SeedingPerformanceWriter::Config seedPerfCfg;
//   seedPerfCfg.inputSeeds = seedingCfg.outputSeeds;
//   seedPerfCfg.inputParticles = inputParticles;
//   seedPerfCfg.inputMeasurementParticlesMap =
//       hitSmearingCfg.outputMeasurementParticlesMap;
//   seedPerfCfg.outputDir = outputDir;
//   seedPerfCfg.outputFilename = "performance_seeding_hists.root";
//   sequencer.addWriter(
//       std::make_shared<SeedingPerformanceWriter>(seedPerfCfg, logLevel));
// 
//   // The track parameters estimation writer
//   RootTrackParameterWriter::Config trackParamsWriterCfg;
//   trackParamsWriterCfg.inputSeeds = seedingCfg.outputSeeds;
//   trackParamsWriterCfg.inputTrackParameters =
//       paramsEstimationCfg.outputTrackParameters;
//   trackParamsWriterCfg.inputTrackParametersSeedMap =
//       paramsEstimationCfg.outputTrackParametersSeedMap;
//   trackParamsWriterCfg.inputParticles = particleReader.outputParticles;
//   trackParamsWriterCfg.inputSimHits = simHitReaderCfg.outputSimHits;
//   trackParamsWriterCfg.inputMeasurementParticlesMap =
//       hitSmearingCfg.outputMeasurementParticlesMap;
//   trackParamsWriterCfg.inputMeasurementSimHitsMap =
//       hitSmearingCfg.outputMeasurementSimHitsMap;
//   trackParamsWriterCfg.outputDir = outputDir;
//   trackParamsWriterCfg.outputFilename = "estimatedparams.root";
//   trackParamsWriterCfg.outputTreename = "estimatedparams";
//   sequencer.addWriter(std::make_shared<RootTrackParameterWriter>(
//       trackParamsWriterCfg, logLevel));
// 
//   return sequencer.run();
}
