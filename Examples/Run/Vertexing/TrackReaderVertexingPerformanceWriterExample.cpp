// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file Find vertices using truth particle information as input
///
/// Reads truth particles from TrackMl files and use the truth information
/// to generate smeared track parameters. Use this pseudo-reconstructed
/// tracks as the input to the vertex finder.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsExamples/Io/Root/RootTrajectorySummaryReader.hpp"
#include "ActsExamples/Io/Root/RootVertexPerformanceWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/ParticleSelectorOptions.hpp"
#include "ActsExamples/Options/VertexingOptions.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"
#include "ActsExamples/TruthTracking/ParticleSelector.hpp"
#include "ActsExamples/TruthTracking/TrackSelector.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithm.hpp"

#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addVertexingOptions(desc);
  Options::addInputOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  Options::addParticleSelectorOptions(desc);
  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // basic setup
  auto logLevel = Options::readLogLevel(vars);
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));
  Sequencer sequencer(Options::readSequencerConfig(vars));

  auto inputDir = vars["input-dir"].as<std::string>();
  auto outputDir =
      ensureWritableDirectory(vars["output-dir"].as<std::string>());

  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vars);

  RootParticleReader::Config particleReaderConfig;
  particleReaderConfig.particleCollection = "allTruthParticles";
  particleReaderConfig.filePath = inputDir + "/particles_initial.root";
  particleReaderConfig.treeName = "particles";
  sequencer.addReader(std::make_shared<RootParticleReader>(
      particleReaderConfig, Acts::Logging::INFO));

  // add additional particle selection
  auto select = Options::readParticleSelectorConfig(vars);
  select.inputParticles = particleReaderConfig.particleCollection;
  select.outputParticles = "detectorAcceptanceSelectedTruthParticles";
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::ParticleSelector>(select, logLevel));

  RootTrajectorySummaryReader::Config trackSummaryReader;
  trackSummaryReader.outputTracks = "fittedTrackParameters";
  trackSummaryReader.outputParticles = "associatedTruthParticles";
  trackSummaryReader.filePath = inputDir + "/tracksummary_fitter.root";
  sequencer.addReader(std::make_shared<RootTrajectorySummaryReader>(
      trackSummaryReader, logLevel));

  // Apply some primary vertexing selection cuts
  TrackSelector::Config trackSelectorConfig;
  trackSelectorConfig.inputTrackParameters = trackSummaryReader.outputTracks;
  trackSelectorConfig.outputTrackParameters = "selectedTracks";
  trackSelectorConfig.removeNeutral = true;
  trackSelectorConfig.absEtaMax = vars["vertexing-eta-max"].as<double>();
  trackSelectorConfig.loc0Max = vars["vertexing-rho-max"].as<double>() * 1_mm;
  trackSelectorConfig.ptMin = vars["vertexing-pt-min"].as<double>() * 1_MeV;
  sequencer.addAlgorithm(
      std::make_shared<TrackSelector>(trackSelectorConfig, logLevel));

  // find vertices
  AdaptiveMultiVertexFinderAlgorithm::Config findVertices;
  findVertices.bField = magneticField;
  findVertices.inputTrackParameters = trackSelectorConfig.outputTrackParameters;
  findVertices.outputProtoVertices = "fittedProtoVertices";
  findVertices.outputVertices = "fittedVertices";
  findVertices.outputTime = "recoTimeMS";
  sequencer.addAlgorithm(std::make_shared<AdaptiveMultiVertexFinderAlgorithm>(
      findVertices, logLevel));

  // write track parameters from fitting
  RootVertexPerformanceWriter::Config vertexWriterConfig;
  vertexWriterConfig.inputAllTruthParticles =
      particleReaderConfig.particleCollection;
  vertexWriterConfig.inputSelectedTruthParticles = select.outputParticles;
  vertexWriterConfig.inputAssociatedTruthParticles =
      trackSummaryReader.outputParticles;
  vertexWriterConfig.inputTrackParameters = trackSummaryReader.outputTracks;
  vertexWriterConfig.inputVertices = findVertices.outputVertices;
  vertexWriterConfig.inputTime = findVertices.outputTime;
  vertexWriterConfig.filePath = outputDir + "/vertexperformance_AMVF.root";
  vertexWriterConfig.treeName = "amvf";
  sequencer.addWriter(std::make_shared<RootVertexPerformanceWriter>(
      vertexWriterConfig, logLevel));

  return sequencer.run();
}
