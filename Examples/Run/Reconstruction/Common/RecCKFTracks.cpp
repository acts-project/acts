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
#include "ActsExamples/Digitization/HitSmearing.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryParametersWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackFindingOptions.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Definitions/Units.hpp>
#include <Acts/MagneticField/BFieldProvider.hpp>

#include <memory>

#include <boost/filesystem.hpp>

#include "RecInput.hpp"

using namespace Acts::UnitLiterals;
using namespace ActsExamples;
using namespace boost::filesystem;
using namespace std::placeholders;

void addRecCKFOptions(ActsExamples::Options::Description& desc) {
  using namespace ActsExamples;
  using boost::program_options::bool_switch;

  auto opt = desc.add_options();
  opt("ckf-truth-seeds", bool_switch(), "Use truth seeds for steering CKF");
}

int runRecCKFTracks(int argc, char* argv[],
                    std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  detector->addOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addTrackFindingOptions(desc);
  addRecCKFOptions(desc);

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
  bool truthSeeded = vm["ckf-truth-seeds"].template as<bool>();

  // Setup detector geometry
  auto geometry = Geometry::build(vm, *detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vm);

  // Read the sim hits
  auto simHitReaderCfg = setupSimHitReading(vm, sequencer);
  // Read the particles
  auto particleReader = setupParticleReading(vm, sequencer);

  // Run the sim hits smearing
  auto hitSmearingCfg = setupSimHitSmearing(
      vm, sequencer, rnd, trackingGeometry, simHitReaderCfg.outputSimHits);

  // Run the particle selection
  // The pre-selection will select truth particles satisfying provided criteria
  // from all particles read in by particle reader for further processing. It
  // has no impact on the truth hits read-in by the cluster reader.
  TruthSeedSelector::Config particleSelectorCfg;
  particleSelectorCfg.inputParticles = particleReader.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.ptMin = 500_MeV;
  particleSelectorCfg.nHitsMin = 9;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  std::string outputTrackParameters;
  // Run particle smearing or run seed finding
  if (truthSeeded) {
    // Run the particle smearing
    auto particleSmearingCfg =
        setupParticleSmearing(vm, sequencer, rnd, inputParticles);
    outputTrackParameters = particleSmearingCfg.outputTrackParameters;
  } else {
    // Create space points
    SpacePointMaker::Config spCfg;
    spCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
    spCfg.inputMeasurements = hitSmearingCfg.outputMeasurements;
    spCfg.outputSpacePoints = "spacepoints";
    spCfg.trackingGeometry = trackingGeometry;
    spCfg.geometrySelection = {
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

    // Algorithm estimating track parameter from seed
    TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
    paramsEstimationCfg.inputSeeds = seedingCfg.outputSeeds;
    paramsEstimationCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
    paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
    paramsEstimationCfg.outputTrackParametersSeedMap =
        "estimatedparams_seed_map";
    paramsEstimationCfg.trackingGeometry = trackingGeometry;
    paramsEstimationCfg.magneticField = magneticField;
    paramsEstimationCfg.bFieldMin = 0.1_T;
    paramsEstimationCfg.sigmaLoc0 = 25._um;
    paramsEstimationCfg.sigmaLoc1 = 100._um;
    paramsEstimationCfg.sigmaPhi = 0.005_degree;
    paramsEstimationCfg.sigmaTheta = 0.001_degree;
    paramsEstimationCfg.sigmaQOverP = 0.1 / 1._GeV;
    paramsEstimationCfg.sigmaT0 = 1400._s;
    sequencer.addAlgorithm(std::make_shared<TrackParamsEstimationAlgorithm>(
        paramsEstimationCfg, logLevel));

    outputTrackParameters = paramsEstimationCfg.outputTrackParameters;
  }

  // Setup the track finding algorithm with CKF
  // It takes all the source links created from truth hit smearing, seeds from
  // truth particle smearing and source link selection config
  auto trackFindingCfg = Options::readTrackFindingConfig(vm);
  trackFindingCfg.inputMeasurements = hitSmearingCfg.outputMeasurements;
  trackFindingCfg.inputSourceLinks = hitSmearingCfg.outputSourceLinks;
  trackFindingCfg.inputInitialTrackParameters = outputTrackParameters;
  trackFindingCfg.outputTrajectories = "trajectories";
  trackFindingCfg.findTracks = TrackFindingAlgorithm::makeTrackFinderFunction(
      trackingGeometry, magneticField);
  sequencer.addAlgorithm(
      std::make_shared<TrackFindingAlgorithm>(trackFindingCfg, logLevel));

  // write track states from CKF
  RootTrajectoryStatesWriter::Config trackStatesWriter;
  trackStatesWriter.inputTrajectories = trackFindingCfg.outputTrajectories;
  // @note The full particles collection is used here to avoid lots of warnings
  // since the unselected CKF track might have a majority particle not in the
  // filtered particle collection. Thsi could be avoided when a seperate track
  // selection algorithm is used.
  trackStatesWriter.inputParticles = particleReader.outputParticles;
  trackStatesWriter.inputSimHits = simHitReaderCfg.outputSimHits;
  trackStatesWriter.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  trackStatesWriter.inputMeasurementSimHitsMap =
      hitSmearingCfg.outputMeasurementSimHitsMap;
  trackStatesWriter.outputDir = outputDir;
  trackStatesWriter.outputFilename = "trackstates_ckf.root";
  trackStatesWriter.outputTreename = "trackstates_ckf";
  sequencer.addWriter(std::make_shared<RootTrajectoryStatesWriter>(
      trackStatesWriter, logLevel));

  // write track parameters from CKF
  RootTrajectoryParametersWriter::Config trackParamsWriter;
  trackParamsWriter.inputTrajectories = trackFindingCfg.outputTrajectories;
  // @note The full particles collection is used here to avoid lots of warnings
  // since the unselected CKF track might have a majority particle not in the
  // filtered particle collection. Thsi could be avoided when a seperate track
  // selection algorithm is used.
  trackParamsWriter.inputParticles = particleReader.outputParticles;
  trackParamsWriter.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  trackParamsWriter.outputDir = outputDir;
  trackParamsWriter.outputFilename = "trackparams_ckf.root";
  trackParamsWriter.outputTreename = "trackparams_ckf";
  sequencer.addWriter(std::make_shared<RootTrajectoryParametersWriter>(
      trackParamsWriter, logLevel));

  // Write CKF performance data
  CKFPerformanceWriter::Config perfWriterCfg;
  perfWriterCfg.inputParticles = inputParticles;
  perfWriterCfg.inputTrajectories = trackFindingCfg.outputTrajectories;
  perfWriterCfg.inputMeasurementParticlesMap =
      hitSmearingCfg.outputMeasurementParticlesMap;
  // The bottom seed on a pixel detector 'eats' two measurements
  perfWriterCfg.nMeasurementsMin = particleSelectorCfg.nHitsMin - 2;
  perfWriterCfg.outputDir = outputDir;
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

  return sequencer.run();
}
