// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "DetectorAlignment.hpp"

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Json/JsonGeometryList.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectorySummaryWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/CsvOptionsReader.hpp"
#include "ActsExamples/Options/DigitizationOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/ParticleSmearingOptions.hpp"
#include "ActsExamples/Options/TrackFittingOptions.hpp"
#include "ActsExamples/Options/TruthSeedSelectorOptions.hpp"
#include "ActsExamples/Reconstruction/ReconstructionBase.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFitterFunction.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Utilities/TracksToTrajectories.hpp"

#include <filesystem>
#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;
using namespace std::filesystem;

void addAlignmentOptions(ActsExamples::Options::Description& desc) {
  using boost::program_options::value;
  auto opt = desc.add_options();
  opt("reco-with-misalignment-correction", value<bool>()->default_value(false),
      "Correct for detector misalignment effects.");
  opt("alignment-geo-config-file", value<std::string>()->default_value(""),
      "Json file for alignment geometry elements selection");
}

int runDetectorAlignment(
    int argc, char* argv[],
    const std::shared_ptr<ActsExamples::IBaseDetector>& detector,
    ActsAlignment::AlignedTransformUpdater alignedTransformUpdater,
    const AlignedDetElementGetter& alignedDetElementsGetter) {
  // using boost::program_options::value;

  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addParticleSmearingOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  detector->addOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addFittingOptions(desc);
  Options::addDigitizationOptions(desc);
  Options::addTruthSeedSelectorOptions(desc);
  addAlignmentOptions(desc);

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

  if (vm["fit-directed-navigation"].as<bool>()) {
    throw std::runtime_error(
        "Directed navigation not supported anymore in the examples binaries."
        "Please refer to the RefittingAlgorithm in the python bindings.");
  }

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
  TruthSeedSelector::Config particleSelectorCfg =
      Options::readTruthSeedSelectorConfig(vm);
  particleSelectorCfg.inputParticles = particleReader.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.nHitsMin = 9;
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

  if (vm["reco-with-misalignment-correction"].as<bool>()) {
    // setup the alignment (which will update the aligned transforms of the
    // detector elements)
    AlignmentAlgorithm::Config alignment;
    alignment.inputSourceLinks = digiCfg.outputSourceLinks;
    alignment.inputMeasurements = digiCfg.outputMeasurements;
    alignment.inputProtoTracks = trackFinderCfg.outputProtoTracks;
    alignment.inputInitialTrackParameters =
        particleSmearingCfg.outputTrackParameters;
    alignment.outputAlignmentParameters = "alignment-parameters";
    alignment.alignedTransformUpdater = std::move(alignedTransformUpdater);
    std::string path = vm["alignment-geo-config-file"].as<std::string>();
    if (not path.empty()) {
      alignment.alignedDetElements = alignedDetElementsGetter(
          detector, ActsExamples::readJsonGeometryList(path));
    }

    // The criteria to determine if the iteration has converged.
    alignment.deltaChi2ONdfCutOff = {10, 0.00005};
    alignment.chi2ONdfCutOff = 0.01;
    alignment.maxNumIterations = 60;
    alignment.align = AlignmentAlgorithm::makeAlignmentFunction(
        trackingGeometry, magneticField);
    sequencer.addAlgorithm(
        std::make_shared<AlignmentAlgorithm>(alignment, logLevel));
  }

  // setup the fitter
  TrackFittingAlgorithm::Config fitter;
  fitter.inputMeasurements = digiCfg.outputMeasurements;
  fitter.inputSourceLinks = digiCfg.outputSourceLinks;
  fitter.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  fitter.inputInitialTrackParameters =
      particleSmearingCfg.outputTrackParameters;
  fitter.outputTracks = "tracks";
  fitter.pickTrack = vm["fit-pick-track"].as<int>();
  fitter.fit = ActsExamples::makeKalmanFitterFunction(
      trackingGeometry, magneticField,
      vm["fit-multiple-scattering-correction"].as<bool>(),
      vm["fit-energy-loss-correction"].as<bool>());
  sequencer.addAlgorithm(
      std::make_shared<TrackFittingAlgorithm>(fitter, logLevel));

  TracksToTrajectories::Config tracksToTrajCfg{};
  tracksToTrajCfg.inputTracks = fitter.outputTracks;
  tracksToTrajCfg.outputTrajectories = "trajectories";
  sequencer.addAlgorithm(
      (std::make_shared<TracksToTrajectories>(tracksToTrajCfg, logLevel)));

  // write track states from fitting
  RootTrajectoryStatesWriter::Config trackStatesWriter;
  trackStatesWriter.inputTrajectories = tracksToTrajCfg.outputTrajectories;
  trackStatesWriter.inputParticles = inputParticles;
  trackStatesWriter.inputSimHits = simHitReaderCfg.outputSimHits;
  trackStatesWriter.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  trackStatesWriter.inputMeasurementSimHitsMap =
      digiCfg.outputMeasurementSimHitsMap;
  trackStatesWriter.filePath = outputDir + "/trackstates_fitter.root";
  sequencer.addWriter(std::make_shared<RootTrajectoryStatesWriter>(
      trackStatesWriter, logLevel));

  // write track summary from CKF
  RootTrajectorySummaryWriter::Config trackSummaryWriter;
  trackSummaryWriter.inputTrajectories = tracksToTrajCfg.outputTrajectories;
  trackSummaryWriter.inputParticles = inputParticles;
  trackSummaryWriter.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  trackSummaryWriter.filePath = outputDir + "/tracksummary_fitter.root";
  sequencer.addWriter(std::make_shared<RootTrajectorySummaryWriter>(
      trackSummaryWriter, logLevel));

  // Write CKF performance data
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
  perfFitter.inputTrajectories = tracksToTrajCfg.outputTrajectories;
  perfFitter.inputParticles = inputParticles;
  perfFitter.inputMeasurementParticlesMap =
      digiCfg.outputMeasurementParticlesMap;
  perfFitter.filePath = outputDir + "/performance_track_fitter.root";
  sequencer.addWriter(
      std::make_shared<TrackFitterPerformanceWriter>(perfFitter, logLevel));

  return sequencer.run();
}
