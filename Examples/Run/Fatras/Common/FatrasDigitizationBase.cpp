// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "FatrasDigitizationBase.hpp"

#include <boost/program_options.hpp>

#include "ACTFW/Digitization/DigitizationAlgorithm.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterWriter.hpp"
#include "ACTFW/Io/Root/RootPlanarClusterWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleStepper.hpp"

void FW::setupDigitization(
    FW::Options::Variables& vars, FW::Sequencer& sequencer,
    std::shared_ptr<const FW::RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry) {
  // Read the standard options
  auto logLevel = FW::Options::readLogLevel(vars);

  // Configure the digitizer
  FW::DigitizationAlgorithm::Config digi;
  digi.inputSimulatedHits = "hits";
  digi.outputClusters = "clusters";
  digi.planarModuleStepper = std::make_shared<Acts::PlanarModuleStepper>(
      Acts::getDefaultLogger("PlanarModuleStepper", logLevel));
  digi.randomNumbers = randomNumbers;
  digi.trackingGeometry = trackingGeometry;
  sequencer.addAlgorithm(
      std::make_shared<FW::DigitizationAlgorithm>(digi, logLevel));

  // Output directory
  std::string outputDir = vars["output-dir"].template as<std::string>();

  // Write digitisation output as Csv files
  if (vars["output-csv"].template as<bool>()) {
    // clusters as root
    FW::CsvPlanarClusterWriter::Config clusterWriterCsv;
    clusterWriterCsv.inputClusters = digi.outputClusters;
    clusterWriterCsv.inputSimulatedHits = digi.inputSimulatedHits;
    clusterWriterCsv.outputDir = outputDir;
    sequencer.addWriter(std::make_shared<FW::CsvPlanarClusterWriter>(
        clusterWriterCsv, logLevel));
  }

  // Write digitsation output as ROOT files
  if (vars["output-root"].template as<bool>()) {
    // clusters as root
    FW::RootPlanarClusterWriter::Config clusterWriterRoot;
    clusterWriterRoot.inputClusters = digi.outputClusters;
    clusterWriterRoot.inputSimulatedHits = digi.inputSimulatedHits;
    clusterWriterRoot.filePath =
        FW::joinPaths(outputDir, digi.outputClusters + ".root");
    sequencer.addWriter(std::make_shared<FW::RootPlanarClusterWriter>(
        clusterWriterRoot, logLevel));
  }
}
