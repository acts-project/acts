// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Digitization/PlanarModuleStepper.hpp"
#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "ActsExamples/Digitization/PlanarSteppingAlgorithm.hpp"
#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterWriter.hpp"
#include "ActsExamples/Io/Root/RootDigitizationWriter.hpp"
#include "ActsExamples/Io/Root/RootPlanarClusterWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <boost/program_options.hpp>

#include "FatrasInternal.hpp"

void addDigitizationOptions(ActsExamples::Options::Description& desc) {
  ActsExamples::Options::addDigitizationOptions(desc);
}

void setupDigitization(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry) {
  using namespace ActsExamples;

  // Read the standard options
  auto logLevel = Options::readLogLevel(vars);
  auto outputDir = vars["output-dir"].template as<std::string>();

  if (vars["digi-smearing"].as<bool>()) {
    SmearingAlgorithm::Config smearCfg = Options::readSmearingConfig(vars);
    smearCfg.inputSimulatedHits = kFatrasCollectionHits;
    smearCfg.outputMeasurements = "measurements";
    smearCfg.trackingGeometry = trackingGeometry;
    smearCfg.randomNumbers = randomNumbers;
    sequencer.addAlgorithm(
        std::make_shared<SmearingAlgorithm>(smearCfg, logLevel));

    // Write digitsation output as ROOT files
    if (vars["output-root"].template as<bool>()) {
      // clusters as root
      RootDigitizationWriter::Config smearWriterRoot;
      smearWriterRoot.inputMeasurements = smearCfg.outputMeasurements;
      smearWriterRoot.inputSimulatedHits = smearCfg.inputSimulatedHits;
      smearWriterRoot.filePath =
          joinPaths(outputDir, smearCfg.outputMeasurements + ".root");
      smearWriterRoot.smearers = smearCfg.smearers;
      sequencer.addWriter(
          std::make_shared<RootDigitizationWriter>(smearWriterRoot, logLevel));
    }

  } else if (vars["digi-geometric-3d"].as<bool>()) {
    // Configure the digitizer
    PlanarSteppingAlgorithm::Config digi;
    digi.inputSimulatedHits = "hits";
    digi.outputClusters = "clusters";
    digi.planarModuleStepper = std::make_shared<Acts::PlanarModuleStepper>(
        Acts::getDefaultLogger("PlanarModuleStepper", logLevel));
    digi.randomNumbers = randomNumbers;
    digi.trackingGeometry = trackingGeometry;
    sequencer.addAlgorithm(
        std::make_shared<PlanarSteppingAlgorithm>(digi, logLevel));

    // Write digitisation output as Csv files
    if (vars["output-csv"].template as<bool>()) {
      // clusters as root
      CsvPlanarClusterWriter::Config clusterWriterCsv;
      clusterWriterCsv.inputClusters = digi.outputClusters;
      clusterWriterCsv.inputSimulatedHits = digi.inputSimulatedHits;
      clusterWriterCsv.outputDir = outputDir;
      sequencer.addWriter(
          std::make_shared<CsvPlanarClusterWriter>(clusterWriterCsv, logLevel));
    }

    // Write digitsation output as ROOT files
    if (vars["output-root"].template as<bool>()) {
      // clusters as root
      RootPlanarClusterWriter::Config clusterWriterRoot;
      clusterWriterRoot.inputClusters = digi.outputClusters;
      clusterWriterRoot.inputSimulatedHits = digi.inputSimulatedHits;
      clusterWriterRoot.filePath =
          joinPaths(outputDir, digi.outputClusters + ".root");
      sequencer.addWriter(std::make_shared<RootPlanarClusterWriter>(
          clusterWriterRoot, logLevel));
    }
  }
}
