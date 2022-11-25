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
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackParameterWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/CsvOptionsReader.hpp"
#include "ActsExamples/Options/DigitizationOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/SpacePointMakerOptions.hpp"
#include "ActsExamples/Reconstruction/ReconstructionBase.hpp"
#include "ActsExamples/TrackFinding/DefaultHoughFunctions.hpp"
#include "ActsExamples/TrackFinding/HoughTransformSeeder.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

#include <boost/program_options.hpp>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int runHoughExample(
    int argc, char* argv[],
    const std::shared_ptr<ActsExamples::IBaseDetector>& detector) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  Options::addInputOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addSpacePointMakerOptions(desc);
  Options::addDigitizationOptions(desc);
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
  for (const auto& cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Setup the magnetic field
  Options::setupMagneticFieldServices(vm, sequencer);
  auto magneticField = Options::readMagneticField(vm);

  // Read the sim hits
  auto simHitReaderCfg = setupSimHitReading(vm, sequencer);
  // Read the particles
  auto particleReader = setupParticleReading(vm, sequencer);

  // Run the sim hits smearing
  auto digiCfg = setupDigitization(vm, sequencer, rnd, tGeometry,
                                   simHitReaderCfg.outputSimHits);

  // Run the particle selection
  // The pre-selection will select truth particles satisfying provided criteria
  // from all particles read in by particle reader for further processing. It
  // has no impact on the truth hits read-in by the cluster reader.
  TruthSeedSelector::Config particleSelectorCfg;
  particleSelectorCfg.inputParticles = particleReader.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.ptMin = 1_GeV;
  particleSelectorCfg.etaMax = 2.5;
  particleSelectorCfg.etaMin = -2.5;
  particleSelectorCfg.nHitsMin = 9;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  // Create space points
  SpacePointMaker::Config spCfg = Options::readSpacePointMakerConfig(vm);
  spCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  spCfg.inputMeasurements = digiCfg.outputMeasurements;
  spCfg.outputSpacePoints = "spacepoints";
  spCfg.trackingGeometry = tGeometry;
  sequencer.addAlgorithm(std::make_shared<SpacePointMaker>(spCfg, logLevel));

  // Hough algorithm
  ActsExamples::HoughTransformSeeder::Config houghCfg;
  houghCfg.inputSpacePoints = {
      spCfg.outputSpacePoints,
  };
  houghCfg.outputSeeds = "seeds";
  houghCfg.outputProtoTracks = "prototracks";

  houghCfg.inputMeasurements = digiCfg.outputMeasurements;
  houghCfg.trackingGeometry = tGeometry;
  houghCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  houghCfg.geometrySelection = {Acts::GeometryIdentifier().setVolume(0)};

  houghCfg.subRegions = {0, 1};

  houghCfg.xMin = -0.05;  // minphi
  houghCfg.xMax = 0.25;   // maxphi
  houghCfg.yMin = -1.25;  // min q/pt, -1/0.8 GeV
  houghCfg.yMax = 1.25;   // max q/pt, +1/0.8 GeV

  houghCfg.houghHistSize_x = 216;
  houghCfg.houghHistSize_y = 216;  // i.e. number of bins in q/pT

  houghCfg.hitExtend_x = {
      2, 1, 0, 0, 0,
      0, 0, 0, 0, 0};  // Hit lines will fill extra bins in x by this amount on
                       // each side, size == nLayers
  houghCfg.nLayers = 10;

  houghCfg.threshold = {9};  // Minimum point value post-convolution to accept
                             // as a road (inclusive)

  houghCfg.localMaxWindowSize =
      0;  // Only create roads from a local maximum, requires traceHits
  houghCfg.kA = 0.0003;  // Assume B = 2T constant.

  houghCfg.fieldCorrector
      .connect<&ActsExamples::DefaultHoughFunctions::fieldCorrectionDefault>();
  houghCfg.layerIDFinder
      .connect<&ActsExamples::DefaultHoughFunctions::findLayerIDDefault>();
  houghCfg.sliceTester
      .connect<&ActsExamples::DefaultHoughFunctions::inSliceDefault>();

  sequencer.addAlgorithm(
      std::make_shared<HoughTransformSeeder>(houghCfg, logLevel));

  // Algorithm estimating track parameter from seed
  TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
  paramsEstimationCfg.inputProtoTracks = houghCfg.outputProtoTracks;
  paramsEstimationCfg.inputSpacePoints = {
      spCfg.outputSpacePoints,
  };
  paramsEstimationCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
  paramsEstimationCfg.outputProtoTracks = "prototracks_estimated";
  paramsEstimationCfg.trackingGeometry = tGeometry;
  paramsEstimationCfg.magneticField = magneticField;
  sequencer.addAlgorithm(std::make_shared<TrackParamsEstimationAlgorithm>(
      paramsEstimationCfg, logLevel));  /// Object 'prototracks' does not exists

  // Seeding performance Writers
  TrackFinderPerformanceWriter::Config tfPerfCfg;
  tfPerfCfg.inputProtoTracks = houghCfg.outputProtoTracks;
  tfPerfCfg.inputParticles = inputParticles;
  tfPerfCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  tfPerfCfg.filePath = outputDir + "/performance_seeding_trees.root";
  sequencer.addWriter(
      std::make_shared<TrackFinderPerformanceWriter>(tfPerfCfg, logLevel));

  SeedingPerformanceWriter::Config seedPerfCfg;
  seedPerfCfg.inputProtoTracks = houghCfg.outputProtoTracks;
  seedPerfCfg.inputParticles = inputParticles;
  seedPerfCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  seedPerfCfg.filePath = outputDir + "/performance_seeding_hists.root";
  sequencer.addWriter(
      std::make_shared<SeedingPerformanceWriter>(seedPerfCfg, logLevel));

  // The track parameters estimation writer
  RootTrackParameterWriter::Config trackParamsWriterCfg;
  trackParamsWriterCfg.inputTrackParameters =
      paramsEstimationCfg.outputTrackParameters;
  trackParamsWriterCfg.inputProtoTracks = paramsEstimationCfg.outputProtoTracks;
  trackParamsWriterCfg.inputParticles = particleReader.outputParticles;
  trackParamsWriterCfg.inputSimHits = simHitReaderCfg.outputSimHits;
  trackParamsWriterCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  trackParamsWriterCfg.inputMeasurementSimHitsMap =
      digiCfg.outputMeasurementSimHitsMap;
  trackParamsWriterCfg.filePath = outputDir + "/estimatedparams.root";
  trackParamsWriterCfg.treeName = "estimatedparams";
  sequencer.addWriter(std::make_shared<RootTrackParameterWriter>(
      trackParamsWriterCfg, logLevel));

  return sequencer.run();
}
