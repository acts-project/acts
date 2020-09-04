// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameters.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Io/Root/RootVertexAndTracksReader.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/Vertexing/IterativeVertexFinderAlgorithm.hpp"

#include <memory>

#include <boost/program_options.hpp>

using namespace ActsExamples;

/// Main executable
///
/// @param argc The argument count
/// @param argv The argument list
int main(int argc, char* argv[]) {
  using namespace boost::program_options;
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addOutputOptions(desc);
  desc.add_options()("input", value<std::string>()->default_value(""),
                     "Input root file to read.");
  auto vm = Options::parse(desc, argc, argv);

  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto logLevel = Options::readLogLevel(vm);

  // Set file to read data from
  std::string fileString = vm["input"].template as<std::string>();

  if (fileString.empty()) {
    std::cout << "Error: Input file not set." << std::endl;
    return EXIT_FAILURE;
  }

  RootVertexAndTracksReader::Config vtxAndTracksReaderCfg;
  vtxAndTracksReaderCfg.fileList.push_back(fileString);

  // Set magnetic field
  Acts::Vector3D bField(0., 0., 2. * Acts::UnitConstants::T);

  // Add the finding algorithm
  ActsExamples::IterativeVertexFinderAlgorithm::Config vertexFindingCfg;
  vertexFindingCfg.trackCollection = vtxAndTracksReaderCfg.outputCollection;
  vertexFindingCfg.bField = bField;

  Sequencer::Config sequencerCfg = Options::readSequencerConfig(vm);
  Sequencer sequencer(sequencerCfg);

  sequencer.addReader(std::make_shared<RootVertexAndTracksReader>(
      vtxAndTracksReaderCfg, logLevel));

  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::IterativeVertexFinderAlgorithm>(
          vertexFindingCfg, logLevel));

  return sequencer.run();
}
