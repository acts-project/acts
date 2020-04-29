// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/EventData/TrackParameters.hpp>
#include <boost/program_options.hpp>
#include <memory>

#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Generators/MultiplicityGenerators.hpp"
#include "ACTFW/Generators/Pythia8ProcessGenerator.hpp"
#include "ACTFW/Generators/VertexGenerators.hpp"
#include "ACTFW/Io/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Io/Root/RootParticleWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Options/Pythia8Options.hpp"
#include "ACTFW/TruthTracking/TrackSelector.hpp"
#include "ACTFW/TruthTracking/TruthVerticesToTracks.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "ACTFW/Vertexing/VertexFitAlgorithm.hpp"

using namespace FW;

/// Main vertex fitter example executable
///
/// @param argc The argument count
/// @param argv The argument list
int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addPythia8Options(desc);
  Options::addOutputOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // basic setup
  auto logLevel = Options::readLogLevel(vm);
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vm));
  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Set up event generator
  EventGenerator::Config evgen = Options::readPythia8Options(vm, logLevel);
  evgen.output = "generated_event";
  evgen.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<EventGenerator>(evgen, logLevel));

  // Set up TruthVerticesToTracks converter algorithm
  TruthVerticesToTracksAlgorithm::Config trkConvConfig;
  trkConvConfig.input = evgen.output;
  trkConvConfig.output = "all_tracks";
  trkConvConfig.randomNumberSvc = rnd;
  trkConvConfig.bField = {0_T, 0_T, 2_T};
  sequencer.addAlgorithm(std::make_shared<TruthVerticesToTracksAlgorithm>(
      trkConvConfig, logLevel));

  // Set up track selector
  TrackSelector::Config selectorConfig;
  selectorConfig.input = trkConvConfig.output;
  selectorConfig.output = "selected_tracks";
  selectorConfig.absEtaMax = 2.5;
  selectorConfig.rhoMax = 4_mm;
  selectorConfig.ptMin = 400_MeV;
  selectorConfig.keepNeutral = false;
  sequencer.addAlgorithm(
      std::make_shared<TrackSelector>(selectorConfig, logLevel));

  // Add the fit algorithm with Billoir fitter
  FWE::VertexFitAlgorithm::Config vertexFitCfg;
  vertexFitCfg.trackCollection = selectorConfig.output;
  vertexFitCfg.bField = trkConvConfig.bField;
  sequencer.addAlgorithm(
      std::make_shared<FWE::VertexFitAlgorithm>(vertexFitCfg, logLevel));

  return sequencer.run();
}
