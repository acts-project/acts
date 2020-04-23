// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <boost/program_options.hpp>

#include "ACTFW/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ACTFW/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Io/Root/RootMaterialTrackWriter.hpp"
#include "ACTFW/Io/Root/RootSimHitWriter.hpp"
#include "ACTFW/MaterialMapping/GeantinoRecording.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/DD4hepG4/DD4hepToG4Svc.hpp"
#include "ACTFW/Utilities/Paths.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  // Declare the supported program options.
  // Setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addOutputOptions(desc);
  FW::Options::addDD4hepOptions(desc);

  // Parse the options
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  FW::Sequencer g4sequencer(FW::Options::readSequencerConfig(vm));

  size_t nTracks = 100;
  int randomSeed1 = 536235167;
  int randomSeed2 = 729237523;

  Acts::GeometryContext geoContext;

  // DETECTOR:
  // --------------------------------------------------------------------------------
  // DD4Hep detector definition
  // read the detector config & dd4hep detector
  auto dd4HepDetectorConfig =
      FW::Options::readDD4hepConfig<po::variables_map>(vm);
  auto geometrySvc =
      std::make_shared<FW::DD4hep::DD4hepGeometryService>(dd4HepDetectorConfig);
  std::shared_ptr<const Acts::TrackingGeometry> tGeometry =
      geometrySvc->trackingGeometry(geoContext);

  // DD4Hep to Geant4 conversion
  //
  FW::DD4hepG4::DD4hepToG4Svc::Config dgConfig("DD4hepToG4",
                                               Acts::Logging::INFO);
  dgConfig.dd4hepService = geometrySvc;
  auto dd4hepToG4Svc = std::make_shared<FW::DD4hepG4::DD4hepToG4Svc>(dgConfig);

  // --------------------------------------------------------------------------------
  // Geant4 JOB:
  // --------------------------------------------------------------------------------
  // set up the writer for
  // ---------------------------------------------------------------------------------

  // set up the algorithm writing out the material map
  FW::GeantinoRecording::Config g4rConfig;
  g4rConfig.geant4Service = dd4hepToG4Svc;
  g4rConfig.tracksPerEvent = nTracks;
  g4rConfig.seed1 = randomSeed1;
  g4rConfig.seed2 = randomSeed2;
  // create the geant4 algorithm
  auto g4rAlgorithm =
      std::make_shared<FW::GeantinoRecording>(g4rConfig, Acts::Logging::INFO);

  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();
  std::string matCollection = g4rConfig.geantMaterialCollection;
  // std::string trkCollection = g4rConfig.geantTrackStepCollection;

  if (vm["output-root"].template as<bool>()) {
    // Write the propagation steps as ROOT TTree
    FW::RootMaterialTrackWriter::Config matTrackWriterRootConfig;
    matTrackWriterRootConfig.prePostStep = true;
    matTrackWriterRootConfig.recalculateTotals = true;
    matTrackWriterRootConfig.collection = matCollection;
    matTrackWriterRootConfig.filePath =
        FW::joinPaths(outputDir, matCollection + ".root");
    auto matTrackWriterRoot =
        std::make_shared<FW::RootMaterialTrackWriter>(matTrackWriterRootConfig);
    g4sequencer.addWriter(matTrackWriterRoot);

    // // Write track step info as ROOT TTree
    // FW::RootSimHitWriter::Config stepSimHitWriterRootConfig;
    // stepSimHitWriterRootConfig.inputSimulatedHits = trkCollection;
    // stepSimHitWriterRootConfig.filePath =
    //     FW::joinPaths(outputDir, trkCollection + ".root");
    // stepSimHitWriterRootConfig.treeName = "steps";
    // auto stepSimHitWriterRoot = std::make_shared<FW::RootSimHitWriter>(
    //     stepSimHitWriterRootConfig, Acts::Logging::INFO);
    // g4sequencer.addWriter(stepSimHitWriterRoot);
  }

  // Append the algorithm and run
  g4sequencer.addAlgorithm(g4rAlgorithm);
  g4sequencer.run();
}
