// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/program_options.hpp>

#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Io/Root/RootMaterialTrackWriter.hpp"
#include "ACTFW/Io/Root/RootSimHitWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "ActsExamples/Geant4/GeantinoRecording.hpp"
#include "G4VUserDetectorConstruction.hh"

using namespace ActsExamples;
using namespace FW;

/// @brief method to process a geometry
/// @param detector The detector descriptor instance
int runSimulation(const boost::program_options::variables_map& vm,
                  std::unique_ptr<G4VUserDetectorConstruction> g4detector) {
  Sequencer sequencer(Options::readSequencerConfig(vm));
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());

  // Setup the Geant4 algorithm
  GeantinoRecording::Config g4;
  std::string materialTrackCollection = g4.outputMaterialTracks;
  g4.detectorConstruction = std::move(g4detector);
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
    materialTrackWriter.collection = materialTrackCollection;
    materialTrackWriter.filePath =
        joinPaths(outputDir, materialTrackCollection + ".root");
    sequencer.addWriter(std::make_shared<RootMaterialTrackWriter>(
        materialTrackWriter, logLevel));
  }
  return sequencer.run();
}