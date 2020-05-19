// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>

#include "ACTFW/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ACTFW/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Io/Root/RootMaterialTrackWriter.hpp"
#include "ACTFW/Io/Root/RootSimHitWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "ActsExamples/Geant4/GeantinoRecording.hpp"
#include "ActsExamples/Geant4DD4hep/DD4hepDetectorConstruction.hpp"

using namespace ActsExamples;
using namespace FW;

int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addOutputOptions(desc);
  Options::addDD4hepOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());

  // setup the DD4hep detector
  auto dd4hepCfg = Options::readDD4hepConfig<po::variables_map>(vm);
  auto geometrySvc = std::make_shared<DD4hep::DD4hepGeometryService>(dd4hepCfg);

  // setup the Geant4 algorithm
  GeantinoRecording::Config g4;
  g4.detectorConstruction =
      std::make_unique<DD4hepDetectorConstruction>(*geometrySvc->lcdd());
  g4.tracksPerEvent = 100;
  g4.seed1 = 536235167;
  g4.seed2 = 729237523;
  sequencer.addAlgorithm(
      std::make_shared<GeantinoRecording>(std::move(g4), logLevel));

  // setup the output writing
  if (vm["output-root"].template as<bool>()) {
    // Write the propagation steps as ROOT TTree
    RootMaterialTrackWriter::Config materialTrackWriter;
    materialTrackWriter.prePostStep = true;
    materialTrackWriter.recalculateTotals = true;
    materialTrackWriter.collection = g4.outputMaterialTracks;
    materialTrackWriter.filePath =
        joinPaths(outputDir, g4.outputMaterialTracks + ".root");
    sequencer.addWriter(std::make_shared<RootMaterialTrackWriter>(
        materialTrackWriter, logLevel));
  }

  return sequencer.run();
}
