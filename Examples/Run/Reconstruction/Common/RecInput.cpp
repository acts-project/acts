// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/IBaseDetector.hpp"
#ifdef ACTS_PLUGIN_ONNX
#include "Acts/Plugins/Onnx/MLTrackClassifier.hpp"
#endif
#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryParametersWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackFindingOptions.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingOptions.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"

#include "RecInput.hpp"

ActsExamples::CsvSimHitReader::Config setupSimHitReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer) {
  using namespace ActsExamples;

  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  // Read truth hits from CSV files
  auto simHitReaderCfg = Options::readCsvSimHitReaderConfig(vars);
  simHitReaderCfg.inputStem = "simhits";
  simHitReaderCfg.outputSimHits = "simhits";
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

ActsExamples::Options::DigitizationConfiguration setupDigitization(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> rnd,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    const std::string& inputSimHits) {
  using namespace ActsExamples;

  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  auto digiCfg = Options::configureDigitization(vars);
  // Common options for digitization
  std::visit([&](auto cfg) {
    cfg.inputSimHits = inputSimHits;
    cfg.outputMeasurements = "measurements";
    cfg.outputSourceLinks = "sourcelinks";
    cfg.outputMeasurementParticlesMap = "measurement_particles_map";
    cfg.outputMeasurementSimHitsMap = "measurement_simhits_map";
    cfg.randomNumbers = rnd;
    cfg.trackingGeometry = trackingGeometry;
  }, digiCfg);
  // Specific to Digitization Algorithm
  if (std::holds_alternative<DigitizationAlgorithm::Config>(digiCfg)) {
    auto cfg = std::get<DigitizationAlgorithm::Config>(digiCfg);
    cfg.outputClusters = "clusters";
    sequencer.addAlgorithm(
      std::make_shared<DigitizationAlgorithm>(cfg, logLevel));
  } else {
    // Specific to SmearingAlgorithm
    auto cfg = std::get<SmearingAlgorithm::Config>(digiCfg);
    sequencer.addAlgorithm(
      std::make_shared<SmearingAlgorithm>(cfg, logLevel));
  }

  return digiCfg;
}

ActsExamples::ParticleSmearing::Config setupParticleSmearing(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> rnd,
    const std::string& inputParticles) {
  using namespace ActsExamples;

  // Read some standard options
  auto logLevel = Options::readLogLevel(vars);

  // Create smeared particles states
  ParticleSmearing::Config particleSmearingCfg;
  particleSmearingCfg.inputParticles = inputParticles;
  particleSmearingCfg.outputTrackParameters = "smearedparameters";
  particleSmearingCfg.randomNumbers = rnd;
  // Gaussian sigmas to smear particle parameters
  particleSmearingCfg.sigmaD0 = 20_um;
  particleSmearingCfg.sigmaD0PtA = 30_um;
  particleSmearingCfg.sigmaD0PtB = 0.3 / 1_GeV;
  particleSmearingCfg.sigmaZ0 = 20_um;
  particleSmearingCfg.sigmaZ0PtA = 30_um;
  particleSmearingCfg.sigmaZ0PtB = 0.3 / 1_GeV;
  particleSmearingCfg.sigmaPhi = 1_degree;
  particleSmearingCfg.sigmaTheta = 1_degree;
  particleSmearingCfg.sigmaPRel = 0.01;
  particleSmearingCfg.sigmaT0 = 1_ns;
  sequencer.addAlgorithm(
      std::make_shared<ParticleSmearing>(particleSmearingCfg, logLevel));

  return particleSmearingCfg;
}
