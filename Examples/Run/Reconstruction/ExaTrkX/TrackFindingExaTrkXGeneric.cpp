// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"  // for evaluating performance
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"  // to read digi config
#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFindingExaTrkX/TrackFindingAlgorithmExaTrkX.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingOptions.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"  // for evaluating performance
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <functional>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

ActsExamples::CsvSimHitReader::Config setupSimHitReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer) {
  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  // Read truth hits from CSV files
  auto simHitReaderCfg = Options::readCsvSimHitReaderConfig(vars);
  simHitReaderCfg.inputStem = "hits";
  simHitReaderCfg.outputSimHits = "hits";
  sequencer.addReader(
      std::make_shared<CsvSimHitReader>(simHitReaderCfg, logLevel));

  return simHitReaderCfg;
}

ActsExamples::CsvParticleReader::Config setupParticleReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer) {
  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  // Read particles (initial states) from CSV files
  auto particleReader = Options::readCsvParticleReaderConfig(vars);
  particleReader.inputStem = "particles_initial";
  particleReader.outputParticles = "particles_initial";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(particleReader, logLevel));

  return particleReader;
}

int main(int argc, char** argv) {
  std::unique_ptr<const Acts::Logger> local_logger;

  // // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::Csv | OutputFormat::Root);
  Options::addInputOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addDigitizationOptions(desc);

  // Add an option for the onnx models
  desc.add_options()("onnx-dir", boost::program_options::value<std::string>(),
                     "The directory where the ONNX models are");

  // Add specific options for this geometry
  // <TODO> make it as an argument.
  auto detector = std::make_shared<GenericDetector>();
  detector->addOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  bool dirNav = false;
  // auto dirNav = vm["directed-navigation"].as<bool>(); // should be false.

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());

  // The geometry, material and decoration
  // build the detector
  auto geometry = Geometry::build(vm, *detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  auto randomNumbers =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vm));

  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Setup the magnetic field
  Options::setupMagneticFieldServices(vm, sequencer);
  auto magneticField = Options::readMagneticField(vm);

  // Read the sim hits
  auto simHitReaderCfg = setupSimHitReading(vm, sequencer);
  // Read the particles
  auto particleReader = setupParticleReading(vm, sequencer);

  // Digitization: SimHits -> Clusters / Measurements
  auto digiCfg = DigitizationConfig(
      vm, readDigiConfigFromJson(vm["digi-config-file"].as<std::string>()));
  digiCfg.inputSimHits = simHitReaderCfg.outputSimHits;
  digiCfg.trackingGeometry = tGeometry;
  digiCfg.randomNumbers = randomNumbers;
  sequencer.addAlgorithm(
      std::make_shared<DigitizationAlgorithm>(digiCfg, logLevel));

  // Select particles in a pre-defined fiducial region for performance study
  TruthSeedSelector::Config particleSelectorCfg;
  particleSelectorCfg.inputParticles = particleReader.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.ptMin = 100_MeV;
  particleSelectorCfg.nHitsMin = 5;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  // measurements:  digiCfg.outputMeasurements
  // clusters: digiCfg.outputClusters
  // Now measurements --> SpacePoints
  SpacePointMaker::Config spCfg;
  spCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  spCfg.inputMeasurements = digiCfg.outputMeasurements;
  spCfg.outputSpacePoints = "spacepoints";
  spCfg.trackingGeometry = tGeometry;
  spCfg.geometrySelection = {Acts::GeometryIdentifier().setVolume(7),
                             Acts::GeometryIdentifier().setVolume(8),
                             Acts::GeometryIdentifier().setVolume(9),
                             Acts::GeometryIdentifier().setVolume(12),
                             Acts::GeometryIdentifier().setVolume(13),
                             Acts::GeometryIdentifier().setVolume(14),
                             Acts::GeometryIdentifier().setVolume(16),
                             Acts::GeometryIdentifier().setVolume(17),
                             Acts::GeometryIdentifier().setVolume(18)};

  sequencer.addAlgorithm(std::make_shared<SpacePointMaker>(spCfg, logLevel));

  // The ML Algorithm

  // ML-based Tracking alg.
  // make module and function names as options..
  TrackFindingAlgorithmExaTrkX::Config trkFinderCfg;
  trkFinderCfg.inputSpacePoints = spCfg.outputSpacePoints;
  trkFinderCfg.outputProtoTracks = "protoTracks";

  if (vm["onnx-dir"].empty() or
      not boost::filesystem::exists(vm["onnx-dir"].as<std::string>())) {
    throw std::invalid_argument(
        "The directory to the ONNX models is empty or invalid");
  }

  Acts::ExaTrkXTrackFinding::Config cfg;
  cfg.inputMLModuleDir = vm["onnx-dir"].as<std::string>();

  trkFinderCfg.trackFinderML = std::make_shared<Acts::ExaTrkXTrackFinding>(cfg);

  sequencer.addAlgorithm(
      std::make_shared<TrackFindingAlgorithmExaTrkX>(trkFinderCfg, logLevel));

  SurfaceSortingAlgorithm::Config sorterCfg;
  // Setup the surface sorter if running direct navigator
  sorterCfg.inputProtoTracks = trkFinderCfg.outputProtoTracks;
  sorterCfg.inputSimulatedHits = simHitReaderCfg.outputSimHits;
  sorterCfg.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  sorterCfg.outputProtoTracks = "sortedprototracks";
  if (dirNav) {
    sequencer.addAlgorithm(
        std::make_shared<SurfaceSortingAlgorithm>(sorterCfg, logLevel));
  }

  auto inputProtoTracks = trkFinderCfg.outputProtoTracks;
  if (dirNav) {
    inputProtoTracks = sorterCfg.outputProtoTracks;
  }

  // Algorithm estimating track parameter from seed
  TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
  paramsEstimationCfg.inputSeeds =
      "";  // it will use spacepoints and input proto tracks as inputs.
  paramsEstimationCfg.inputProtoTracks = inputProtoTracks;
  paramsEstimationCfg.inputSpacePoints = {
      spCfg.outputSpacePoints,
  };
  paramsEstimationCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
  paramsEstimationCfg.outputProtoTracks = "prototracks_estimated";
  paramsEstimationCfg.trackingGeometry = tGeometry;
  paramsEstimationCfg.magneticField = magneticField;
  paramsEstimationCfg.bFieldMin = 0.1_T;
  paramsEstimationCfg.deltaRMax = 100._mm;
  paramsEstimationCfg.sigmaLoc0 = 25._um;
  paramsEstimationCfg.sigmaLoc1 = 100._um;
  paramsEstimationCfg.sigmaPhi = 0.005_degree;
  paramsEstimationCfg.sigmaTheta = 0.001_degree;
  paramsEstimationCfg.sigmaQOverP = 0.1 / 1._GeV;
  paramsEstimationCfg.sigmaT0 = 1400._s;
  sequencer.addAlgorithm(std::make_shared<TrackParamsEstimationAlgorithm>(
      paramsEstimationCfg, logLevel));

  // Track fitting
  // setup the fitter
  TrackFittingAlgorithm::Config fitter;
  fitter.inputMeasurements = digiCfg.outputMeasurements;
  fitter.inputSourceLinks = digiCfg.outputSourceLinks;
  fitter.inputProtoTracks = trkFinderCfg.outputProtoTracks;
  if (dirNav) {
    fitter.inputProtoTracks = sorterCfg.outputProtoTracks;
  }
  fitter.inputInitialTrackParameters =
      paramsEstimationCfg.outputTrackParameters;
  fitter.outputTrajectories = "trajectories";
  fitter.directNavigation = dirNav;
  fitter.trackingGeometry = tGeometry;
  fitter.dFit = TrackFittingAlgorithm::makeKalmanFitterFunction(magneticField);
  fitter.fit =
      TrackFittingAlgorithm::makeKalmanFitterFunction(tGeometry, magneticField);
  sequencer.addAlgorithm(
      std::make_shared<TrackFittingAlgorithm>(fitter, logLevel));

  // write out performance
  // write track finding/seeding performance
  TrackFinderPerformanceWriter::Config tfPerfCfg;
  tfPerfCfg.inputProtoTracks = trkFinderCfg.outputProtoTracks;
  // using selected particles
  tfPerfCfg.inputParticles = inputParticles;
  tfPerfCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  sequencer.addWriter(
      std::make_shared<TrackFinderPerformanceWriter>(tfPerfCfg, logLevel));

  // Write track finding performance data
  CKFPerformanceWriter::Config perfWriterCfg;
  perfWriterCfg.inputParticles = inputParticles;
  perfWriterCfg.inputTrajectories = fitter.outputTrajectories;
  perfWriterCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  // The bottom seed on a pixel detector 'eats' one or two measurements?
  perfWriterCfg.nMeasurementsMin = particleSelectorCfg.nHitsMin;
  sequencer.addWriter(
      std::make_shared<CKFPerformanceWriter>(perfWriterCfg, logLevel));

  return sequencer.run();
}
