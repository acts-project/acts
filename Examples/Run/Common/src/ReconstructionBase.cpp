// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Reconstruction/ReconstructionBase.hpp"

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/DigitizationOptions.hpp"
#include "ActsExamples/Options/TrackFindingOptions.hpp"
#include "ActsExamples/Options/TrackFittingOptions.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Options.hpp"

ActsExamples::CsvSimHitReader::Config setupSimHitReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer) {
  using namespace ActsExamples;

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
  using namespace ActsExamples;

  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  // Read particles (initial states) and clusters from CSV files
  auto particleReader = Options::readCsvParticleReaderConfig(vars);
  particleReader.inputStem = "particles_initial";
  particleReader.outputParticles = "particles_initial";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(particleReader, logLevel));

  return particleReader;
}

ActsExamples::DigitizationConfig setupDigitization(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> rnd,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    const std::string& inputSimHits) {
  using namespace ActsExamples;

  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  auto digiCfg = ActsExamples::DigitizationConfig(
      vars["digi-merge"].as<bool>(), vars["digi-merge-nsigma"].as<double>(),
      vars["digi-merge-common-corner"].as<bool>(),
      ActsExamples::readDigiConfigFromJson(
          vars["digi-config-file"].as<std::string>()));
  // Common options for digitization
  digiCfg.inputSimHits = inputSimHits;
  digiCfg.randomNumbers = std::move(rnd);
  digiCfg.trackingGeometry = std::move(trackingGeometry);
  sequencer.addAlgorithm(
      std::make_shared<DigitizationAlgorithm>(digiCfg, logLevel));

  if (not vars["dump-digi-config"].as<std::string>().empty()) {
    writeDigiConfigToJson(digiCfg.digitizationConfigs,
                          vars["dump-digi-config"].as<std::string>());
  }

  return digiCfg;
}

ActsExamples::ParticleSmearing::Config setupParticleSmearing(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> rnd,
    const std::string& inputParticles) {
  using namespace ActsExamples;
  using namespace Acts::UnitLiterals;

  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  // Create smeared particles states
  ParticleSmearing::Config particleSmearingCfg =
      Options::readParticleSmearingOptions(vars);
  particleSmearingCfg.inputParticles = inputParticles;
  particleSmearingCfg.outputTrackParameters = "smearedparameters";
  particleSmearingCfg.randomNumbers = std::move(rnd);
  sequencer.addAlgorithm(
      std::make_shared<ParticleSmearing>(particleSmearingCfg, logLevel));

  return particleSmearingCfg;
}

ActsExamples::CsvMeasurementReader::Config setupMeasurementReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer) {
  using namespace ActsExamples;
  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  // Read particles (initial states) from CSV files
  auto measurementsReader = Options::readCsvMeasurementReaderConfig(vars);
  measurementsReader.outputMeasurements = "measurements";
  measurementsReader.outputMeasurementSimHitsMap = "measurements2hits";
  measurementsReader.outputSourceLinks = "source_links";
  measurementsReader.outputClusters = "clusters";
  sequencer.addReader(
      std::make_shared<CsvMeasurementReader>(measurementsReader, logLevel));

  return measurementsReader;
}
