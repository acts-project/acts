// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/IBaseDetector.hpp"
#ifdef ACTS_PLUGIN_ONNX
#include "Acts/Plugins/Onnx/MLTrackClassifier.hpp"
#endif
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvMultiTrajectoryWriter.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectorySummaryWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/CsvOptionsReader.hpp"
#include "ActsExamples/Options/CsvOptionsWriter.hpp"
#include "ActsExamples/Options/DigitizationOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/SpacePointMakerOptions.hpp"
#include "ActsExamples/Options/TrackFindingOptions.hpp"
#include "ActsExamples/Options/TrackFittingOptions.hpp"
#include "ActsExamples/Reconstruction/ReconstructionBase.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/SeedsToPrototracks.hpp"
#include "ActsExamples/Utilities/TracksToTrajectories.hpp"
#include <Acts/Definitions/Units.hpp>

#include <filesystem>
#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;
using namespace std::filesystem;
using namespace std::placeholders;

void addRecCKFOptions(ActsExamples::Options::Description& desc) {
  using namespace ActsExamples;
  using boost::program_options::bool_switch;

  auto opt = desc.add_options();
  opt("ckf-truth-smeared-seeds", bool_switch(),
      "Use track parameters smeared from truth particles for steering CKF");
  opt("ckf-truth-estimated-seeds", bool_switch(),
      "Use track parameters estimated from truth tracks for steering CKF");
}

int runRecCKFTracks(
    int argc, char* argv[],
    const std::shared_ptr<ActsExamples::IBaseDetector>& detector) {
  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc,
                            OutputFormat::Csv | OutputFormat::DirectoryOnly);
  detector->addOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addFittingOptions(desc);
  Options::addTrackFindingOptions(desc);
  addRecCKFOptions(desc);
  Options::addDigitizationOptions(desc);
  Options::addParticleSmearingOptions(desc);
  Options::addSpacePointMakerOptions(desc);
  Options::addCsvWriterOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto inputDir = vm["input-dir"].as<std::string>();
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(
      Options::readRandomNumbersConfig(vm));
  bool truthSmearedSeeded = vm["ckf-truth-smeared-seeds"].template as<bool>();
  bool truthEstimatedSeeded =
      vm["ckf-truth-estimated-seeds"].template as<bool>();

  // Setup detector geometry
  auto geometry = Geometry::build(vm, *detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (const auto& cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vm);

  // Read the sim hits
  auto simHitReaderCfg = setupSimHitReading(vm, sequencer);
  // Read the particles
  auto particleReader = setupParticleReading(vm, sequencer);

  // Run the sim hits smearing
  auto digiCfg = setupDigitization(vm, sequencer, rnd, trackingGeometry,
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
  particleSelectorCfg.ptMin = 500_MeV;
  particleSelectorCfg.nHitsMin = 9;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  // Create starting parameters from either particle smearing or combined seed
  // finding and track parameters estimation
  std::string outputTrackParameters;
  if (truthSmearedSeeded) {
    // Run the particle smearing
    auto particleSmearingCfg =
        setupParticleSmearing(vm, sequencer, rnd, inputParticles);
    outputTrackParameters = particleSmearingCfg.outputTrackParameters;
  } else {
    // Create space points
    SpacePointMaker::Config spCfg = Options::readSpacePointMakerConfig(vm);
    spCfg.inputSourceLinks = digiCfg.outputSourceLinks;
    spCfg.inputMeasurements = digiCfg.outputMeasurements;
    spCfg.outputSpacePoints = "spacepoints";
    spCfg.trackingGeometry = trackingGeometry;
    sequencer.addAlgorithm(std::make_shared<SpacePointMaker>(spCfg, logLevel));

    // Create seeds (i.e. proto tracks) using either truth track finding or seed
    // finding algorithm
    std::string inputProtoTracks = "";
    std::string inputSeeds = "";
    if (truthEstimatedSeeded) {
      // Truth track finding algorithm
      TruthTrackFinder::Config trackFinderCfg;
      trackFinderCfg.inputParticles = inputParticles;
      trackFinderCfg.inputMeasurementParticlesMap =
          digiCfg.outputMeasurementParticlesMap;
      trackFinderCfg.outputProtoTracks = "prototracks";
      sequencer.addAlgorithm(
          std::make_shared<TruthTrackFinder>(trackFinderCfg, logLevel));
      inputProtoTracks = trackFinderCfg.outputProtoTracks;
    } else {
      // Seeding algorithm
      SeedingAlgorithm::Config seedingCfg;
      seedingCfg.inputSpacePoints = {
          spCfg.outputSpacePoints,
      };
      seedingCfg.outputSeeds = "seeds";

      seedingCfg.gridConfig.rMax = 200._mm;
      seedingCfg.seedFinderConfig.rMax = seedingCfg.gridConfig.rMax;

      seedingCfg.seedFilterConfig.deltaRMin = 1_mm;
      seedingCfg.seedFinderConfig.deltaRMin =
          seedingCfg.seedFilterConfig.deltaRMin;

      seedingCfg.gridConfig.deltaRMax = 60._mm;
      seedingCfg.seedFinderConfig.deltaRMax = seedingCfg.gridConfig.deltaRMax;

      seedingCfg.seedFinderConfig.collisionRegionMin = -250_mm;
      seedingCfg.seedFinderConfig.collisionRegionMax = 250._mm;

      seedingCfg.gridConfig.zMin = -2000._mm;
      seedingCfg.gridConfig.zMax = 2000._mm;
      seedingCfg.seedFinderConfig.zMin = seedingCfg.gridConfig.zMin;
      seedingCfg.seedFinderConfig.zMax = seedingCfg.gridConfig.zMax;

      seedingCfg.seedFilterConfig.maxSeedsPerSpM = 1;
      seedingCfg.seedFinderConfig.maxSeedsPerSpM =
          seedingCfg.seedFilterConfig.maxSeedsPerSpM;

      seedingCfg.gridConfig.cotThetaMax = 7.40627;  // 2.7 eta
      seedingCfg.seedFinderConfig.cotThetaMax =
          seedingCfg.gridConfig.cotThetaMax;

      seedingCfg.seedFinderConfig.sigmaScattering = 50;
      seedingCfg.seedFinderConfig.radLengthPerSeed = 0.1;

      seedingCfg.gridConfig.minPt = 500._MeV;
      seedingCfg.seedFinderConfig.minPt = seedingCfg.gridConfig.minPt;

      seedingCfg.gridOptions.bFieldInZ = 1.99724_T;

      seedingCfg.seedFinderOptions.bFieldInZ = seedingCfg.gridOptions.bFieldInZ;
      seedingCfg.seedFinderOptions.beamPos = {0_mm, 0_mm};

      seedingCfg.seedFinderConfig.impactMax = 3._mm;

      sequencer.addAlgorithm(
          std::make_shared<SeedingAlgorithm>(seedingCfg, logLevel));

      SeedsToPrototracks::Config seedsToPrototrackCfg;
      seedsToPrototrackCfg.inputSeeds = seedingCfg.outputSeeds;
      seedsToPrototrackCfg.outputProtoTracks = "prototracks";
      sequencer.addAlgorithm(
          std::make_shared<SeedsToPrototracks>(seedsToPrototrackCfg, logLevel));

      inputProtoTracks = seedsToPrototrackCfg.outputProtoTracks;
      inputSeeds = seedingCfg.outputSeeds;
    }

    // write track finding/seeding performance
    TrackFinderPerformanceWriter::Config tfPerfCfg;
    tfPerfCfg.inputProtoTracks = inputProtoTracks;
    // using selected particles
    tfPerfCfg.inputParticles = inputParticles;
    tfPerfCfg.inputMeasurementParticlesMap =
        digiCfg.outputMeasurementParticlesMap;
    tfPerfCfg.filePath = outputDir + "/performance_seeding_trees.root";
    sequencer.addWriter(
        std::make_shared<TrackFinderPerformanceWriter>(tfPerfCfg, logLevel));

    // Algorithm estimating track parameter from seed
    TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
    paramsEstimationCfg.inputSeeds = inputSeeds;
    paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
    paramsEstimationCfg.trackingGeometry = trackingGeometry;
    paramsEstimationCfg.magneticField = magneticField;
    paramsEstimationCfg.bFieldMin = 0.1_T;
    paramsEstimationCfg.sigmaLoc0 = 25._um;
    paramsEstimationCfg.sigmaLoc1 = 100._um;
    paramsEstimationCfg.sigmaPhi = 0.02_degree;
    paramsEstimationCfg.sigmaTheta = 0.02_degree;
    paramsEstimationCfg.sigmaQOverP = 0.1 / 1._GeV;
    paramsEstimationCfg.sigmaT0 = 1400._s;
    paramsEstimationCfg.initialVarInflation =
        vm["ckf-initial-variance-inflation"].template as<Options::Reals<6>>();

    sequencer.addAlgorithm(std::make_shared<TrackParamsEstimationAlgorithm>(
        paramsEstimationCfg, logLevel));

    outputTrackParameters = paramsEstimationCfg.outputTrackParameters;
  }

  // Setup the track finding algorithm with CKF
  // It takes all the source links created from truth hit smearing, seeds from
  // truth particle smearing and source link selection config
  auto trackFindingCfg = Options::readTrackFindingConfig(vm);
  trackFindingCfg.inputMeasurements = digiCfg.outputMeasurements;
  trackFindingCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  trackFindingCfg.inputInitialTrackParameters = outputTrackParameters;
  trackFindingCfg.outputTracks = "tracks";
  trackFindingCfg.computeSharedHits = true;
  trackFindingCfg.findTracks = TrackFindingAlgorithm::makeTrackFinderFunction(
      trackingGeometry, magneticField);
  sequencer.addAlgorithm(
      std::make_shared<TrackFindingAlgorithm>(trackFindingCfg, logLevel));

  TracksToTrajectories::Config tracksToTrajCfg{};
  tracksToTrajCfg.inputTracks = trackFindingCfg.outputTracks;
  tracksToTrajCfg.outputTrajectories = "trajectories";
  sequencer.addAlgorithm(
      (std::make_shared<TracksToTrajectories>(tracksToTrajCfg, logLevel)));

  // write track states from CKF
  RootTrajectoryStatesWriter::Config trackStatesWriter;
  trackStatesWriter.inputTrajectories = tracksToTrajCfg.outputTrajectories;
  // @note The full particles collection is used here to avoid lots of warnings
  // since the unselected CKF track might have a majority particle not in the
  // filtered particle collection. This could be avoided when a seperate track
  // selection algorithm is used.
  trackStatesWriter.inputParticles = particleReader.outputParticles;
  trackStatesWriter.inputSimHits = simHitReaderCfg.outputSimHits;
  trackStatesWriter.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  trackStatesWriter.inputMeasurementSimHitsMap =
      digiCfg.outputMeasurementSimHitsMap;
  trackStatesWriter.filePath = outputDir + "/trackstates_ckf.root";
  trackStatesWriter.treeName = "trackstates";
  sequencer.addWriter(std::make_shared<RootTrajectoryStatesWriter>(
      trackStatesWriter, logLevel));

  // write track summary from CKF
  RootTrajectorySummaryWriter::Config trackSummaryWriter;
  trackSummaryWriter.inputTrajectories = tracksToTrajCfg.outputTrajectories;
  // @note The full particles collection is used here to avoid lots of warnings
  // since the unselected CKF track might have a majority particle not in the
  // filtered particle collection. This could be avoided when a seperate track
  // selection algorithm is used.
  trackSummaryWriter.inputParticles = particleReader.outputParticles;
  trackSummaryWriter.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  trackSummaryWriter.filePath = outputDir + "/tracksummary_ckf.root";
  trackSummaryWriter.treeName = "tracksummary";
  sequencer.addWriter(std::make_shared<RootTrajectorySummaryWriter>(
      trackSummaryWriter, logLevel));

  // Write CKF performance data
  CKFPerformanceWriter::Config perfWriterCfg;
  perfWriterCfg.inputParticles = inputParticles;
  perfWriterCfg.inputTrajectories = tracksToTrajCfg.outputTrajectories;
  perfWriterCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  // The bottom seed could be the first, second or third hits on the truth track
  perfWriterCfg.nMeasurementsMin = particleSelectorCfg.nHitsMin - 3;
  perfWriterCfg.ptMin = 0.4_GeV;
  perfWriterCfg.filePath = outputDir + "/performance_ckf.root";
#ifdef ACTS_PLUGIN_ONNX
  // Onnx plugin related options
  // Path to default demo ML model for track classification
  path currentFilePath(__FILE__);
  path parentPath = currentFilePath.parent_path();
  path demoModelPath =
      canonical(parentPath / "MLAmbiguityResolutionDemo.onnx").native();
  // Threshold probability for neural network to classify track as duplicate
  double decisionThreshProb = 0.5;
  // Initialize OnnxRuntime plugin
  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "MLTrackClassifier");
  Acts::MLTrackClassifier neuralNetworkClassifier(env, demoModelPath.c_str());
  perfWriterCfg.duplicatedPredictor =
      std::bind(&Acts::MLTrackClassifier::isDuplicate, &neuralNetworkClassifier,
                std::placeholders::_1, decisionThreshProb);
#endif
  sequencer.addWriter(
      std::make_shared<CKFPerformanceWriter>(perfWriterCfg, logLevel));

  if (vm["output-csv"].template as<bool>()) {
    // Write the CKF track as Csv
    CsvMultiTrajectoryWriter::Config trackWriterCsvConfig;
    trackWriterCsvConfig.inputTrajectories = tracksToTrajCfg.outputTrajectories;
    trackWriterCsvConfig.outputDir = outputDir;
    trackWriterCsvConfig.inputMeasurementParticlesMap =
        digiCfg.outputMeasurementParticlesMap;
    sequencer.addWriter(std::make_shared<CsvMultiTrajectoryWriter>(
        trackWriterCsvConfig, logLevel));
  }

  return sequencer.run();
}
