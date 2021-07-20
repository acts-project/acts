// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geant4/G4DetectorConstructionFactory.hpp"
#include "ActsExamples/Geant4/Geant4Options.hpp"
#include "ActsExamples/Geant4/GeantinoRecording.hpp"
#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <boost/program_options.hpp>

#include "G4VUserDetectorConstruction.hh"

/// @brief method to process a geometry
/// @param detector The detector descriptor instance
inline int runGeantinoRecording(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<ActsExamples::G4DetectorConstructionFactory>
        g4DetectorFactory) {
  using namespace ActsExamples;
  Sequencer sequencer(Options::readSequencerConfig(vm));
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());

  // Setup the Geant4 algorithm
  auto g4Config = Options::readGeantinoRecordingConfig(vm);
  auto outputMaterialTracks = g4Config.outputMaterialTracks;
  g4Config.detectorConstructionFactory = std::move(g4DetectorFactory);
  sequencer.addAlgorithm(
      std::make_shared<GeantinoRecording>(std::move(g4Config), logLevel));

  // setup the output writing
  if (vm["output-root"].template as<bool>()) {
    // Write the propagation steps as ROOT TTree
    RootMaterialTrackWriter::Config materialTrackWriter;
    materialTrackWriter.prePostStep = true;
    materialTrackWriter.recalculateTotals = true;
    materialTrackWriter.collection = outputMaterialTracks;
    materialTrackWriter.filePath =
        joinPaths(outputDir, outputMaterialTracks + ".root");
    sequencer.addWriter(std::make_shared<RootMaterialTrackWriter>(
        materialTrackWriter, logLevel));
  }
  return sequencer.run();
}