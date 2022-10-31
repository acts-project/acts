// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"  // TODO: optimize for chi2?
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectorySummaryWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/CsvOptionsReader.hpp"
#include "ActsExamples/Options/DigitizationOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/TrackFittingChi2Options.hpp"
#include "ActsExamples/Options/TrackFittingOptions.hpp"
#include "ActsExamples/Options/TruthSeedSelectorOptions.hpp"
#include "ActsExamples/Reconstruction/ReconstructionBase.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFittingChi2/TrackFittingChi2Algorithm.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int runRecChi2Tracks(int argc, char* argv[],
                     std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  using boost::program_options::value;

  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  detector->addOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addFittingOptions(desc);  // TODO: not all options apply to chi2.
                                     // copy them over to Chi2Options instead ?
  Options::addFittingChi2Options(desc);
  Options::addDigitizationOptions(desc);
  Options::addParticleSmearingOptions(desc);
  Options::addTruthSeedSelectorOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  auto rnd = std::make_shared<const ActsExamples::RandomNumbers>(
      Options::readRandomNumbersConfig(vm));
  // auto dirNav = vm["fit-directed-navigation"].as<bool>(); // TODO: not used
  auto nUpdates = vm["chi2-updates"].as<unsigned int>();

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
  auto digiCfg = setupDigitization(vm, sequencer, rnd, trackingGeometry,
                                   simHitReaderCfg.outputSimHits);
  // Run the particle selection
  // The pre-selection will select truth particles satisfying provided criteria
  // from all particles read in by particle reader for further processing. It
  // has no impact on the truth hits read-in by the cluster reader.
  TruthSeedSelector::Config particleSelectorCfg =
      Options::readTruthSeedSelectorConfig(vm);
  particleSelectorCfg.inputParticles = particleReader.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.nHitsMin = 9;  // 9-layer telescope
  particleSelectorCfg.ptMin = 500._MeV;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  // Run the particle smearing
  auto particleSmearingCfg =
      setupParticleSmearing(vm, sequencer, rnd, inputParticles);

  // The fitter needs the measurements (proto tracks) and initial
  // track states (proto states). The elements in both collections
  // must match and must be created from the same input particles.
  // Create truth tracks
  TruthTrackFinder::Config trackFinderCfg;
  trackFinderCfg.inputParticles = inputParticles;
  trackFinderCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  trackFinderCfg.outputProtoTracks = "prototracks";
  sequencer.addAlgorithm(
      std::make_shared<TruthTrackFinder>(trackFinderCfg, logLevel));

  //   SurfaceSortingAlgorithm::Config sorterCfg;
  //   // Setup the surface sorter if running direct navigator
  //   sorterCfg.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  //   sorterCfg.inputSimHits = simHitReaderCfg.outputSimHits;
  //   sorterCfg.inputMeasurementSimHitsMap =
  //   digiCfg.outputMeasurementSimHitsMap; sorterCfg.outputProtoTracks =
  //   "sortedprototracks"; if (dirNav) {
  //     sequencer.addAlgorithm(
  //         std::make_shared<SurfaceSortingAlgorithm>(sorterCfg, logLevel));
  //   }

  // setup the fitter
  // TODO: wouldn't it be better to pass some chi2 options directly here,
  // instead of constructing the Acts::Chi2FitterOptions object in the
  // TrackFittingAlgorithm.cpp?
  TrackFittingChi2Algorithm::Config fitter;
  fitter.inputMeasurements = digiCfg.outputMeasurements;
  fitter.inputSourceLinks = digiCfg.outputSourceLinks;
  fitter.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  // if (dirNav) {
  //   fitter.inputProtoTracks = sorterCfg.outputProtoTracks;
  // }
  fitter.inputInitialTrackParameters =
      particleSmearingCfg.outputTrackParameters;
  fitter.outputTrajectories = "trajectories";
  // fitter.directNavigation = false; // TODO: the --direct-navigation flag is
  // ignored for now.
  fitter.multipleScattering =
      vm["fit-multiple-scattering-correction"].as<bool>();
  fitter.energyLoss = vm["fit-energy-loss-correction"].as<bool>();

  fitter.multipleScattering = false;  // TODO
  fitter.energyLoss = false;          // TODO
  // ^ not implemented yet. both have a default value of 'true' in the options,
  // so we hard-code them here.

  fitter.pickTrack = vm["fit-pick-track"].as<int>();
  fitter.trackingGeometry = trackingGeometry;
  fitter.nUpdates = nUpdates;
  fitter.fit = TrackFittingChi2Algorithm::makeTrackFitterChi2Function(
      trackingGeometry, magneticField);
  // ^ returns a Chi2TrackFitterFunction
  sequencer.addAlgorithm(
      std::make_shared<TrackFittingChi2Algorithm>(fitter, logLevel));

  // write track states from fitting
  RootTrajectoryStatesWriter::Config trackStatesWriter;
  trackStatesWriter.inputTrajectories = fitter.outputTrajectories;
  trackStatesWriter.inputParticles = inputParticles;
  trackStatesWriter.inputSimHits = simHitReaderCfg.outputSimHits;
  trackStatesWriter.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  trackStatesWriter.inputMeasurementSimHitsMap =
      digiCfg.outputMeasurementSimHitsMap;
  trackStatesWriter.filePath = outputDir + "/trackstates_fitter.root";
  sequencer.addWriter(std::make_shared<RootTrajectoryStatesWriter>(
      trackStatesWriter, logLevel));

  // write track summary
  RootTrajectorySummaryWriter::Config trackSummaryWriter;
  trackSummaryWriter.inputTrajectories = fitter.outputTrajectories;
  trackSummaryWriter.inputParticles = inputParticles;
  trackSummaryWriter.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  trackSummaryWriter.filePath = outputDir + "/tracksummary_fitter.root";
  sequencer.addWriter(std::make_shared<RootTrajectorySummaryWriter>(
      trackSummaryWriter, logLevel));

  // write reconstruction performance data
  TrackFinderPerformanceWriter::Config perfFinder;
  perfFinder.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  perfFinder.inputParticles = inputParticles;
  perfFinder.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  perfFinder.filePath = outputDir + "/performance_track_finder.root";
  sequencer.addWriter(
      std::make_shared<TrackFinderPerformanceWriter>(perfFinder, logLevel));

  TrackFitterPerformanceWriter::Config perfFitter;
  perfFitter.inputTrajectories = fitter.outputTrajectories;
  perfFitter.inputParticles = inputParticles;
  perfFitter.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  perfFitter.filePath = outputDir + "/performance_track_fitter.root";
  sequencer.addWriter(
      std::make_shared<TrackFitterPerformanceWriter>(perfFitter, logLevel));

  return sequencer.run();
}
