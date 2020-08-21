// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "FatrasDigitizationBase.hpp"

#include "Acts/Plugins/Digitization/PlanarModuleStepper.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterWriter.hpp"
#include "ActsExamples/Io/Root/RootPlanarClusterWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <boost/program_options.hpp>

void ActsExamples::setupDigitization(
    ActsExamples::Options::Variables& vars, ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry) {
  // Read the standard options
  auto logLevel = ActsExamples::Options::readLogLevel(vars);

  // Configure the digitizer
  ActsExamples::DigitizationAlgorithm::Config digi;
  digi.inputSimulatedHits = "hits";
  digi.outputClusters = "clusters";
  digi.planarModuleStepper = std::make_shared<Acts::PlanarModuleStepper>(
      Acts::getDefaultLogger("PlanarModuleStepper", logLevel));
  digi.randomNumbers = randomNumbers;
  digi.trackingGeometry = trackingGeometry;
  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::DigitizationAlgorithm>(digi, logLevel));

  // Output directory
  std::string outputDir = vars["output-dir"].template as<std::string>();

  // Write digitisation output as Csv files
  if (vars["output-csv"].template as<bool>()) {
    // clusters as root
    ActsExamples::CsvPlanarClusterWriter::Config clusterWriterCsv;
    clusterWriterCsv.inputClusters = digi.outputClusters;
    clusterWriterCsv.inputSimulatedHits = digi.inputSimulatedHits;
    clusterWriterCsv.outputDir = outputDir;
    sequencer.addWriter(std::make_shared<ActsExamples::CsvPlanarClusterWriter>(
        clusterWriterCsv, logLevel));
  }

  // Write digitsation output as ROOT files
  if (vars["output-root"].template as<bool>()) {
    // clusters as root
    ActsExamples::RootPlanarClusterWriter::Config clusterWriterRoot;
    clusterWriterRoot.inputClusters = digi.outputClusters;
    clusterWriterRoot.inputSimulatedHits = digi.inputSimulatedHits;
    clusterWriterRoot.filePath =
        ActsExamples::joinPaths(outputDir, digi.outputClusters + ".root");
    sequencer.addWriter(std::make_shared<ActsExamples::RootPlanarClusterWriter>(
        clusterWriterRoot, logLevel));
  }
}
