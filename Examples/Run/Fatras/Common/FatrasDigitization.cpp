// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleStepper.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "ActsExamples/Digitization/PlanarSteppingAlgorithm.hpp"
#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterWriter.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Io/Root/RootDigitizationWriter.hpp"
#include "ActsExamples/Io/Root/RootPlanarClusterWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <fstream>

#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

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

  auto cfile = vars["digi-config-file"].as<std::string>();

  if (not cfile.empty()) {
    auto in = std::ifstream(cfile, std::ifstream::in | std::ifstream::binary);
    if (in.good()) {
      // Get the json file
      nlohmann::json djson;
      in >> djson;
      in.close();

      DigiConfigContainer digitizationConfigs =
          DigiConfigConverter("digitization-configuration").fromJson(djson);

      DigitizationAlgorithm::Config digiCfg =
          Options::readDigitizationConfig(vars);
      digiCfg.inputSimHits = kFatrasCollectionHits;
      digiCfg.outputMeasurements = "measurements";
      digiCfg.outputSourceLinks = "sourcelinks";
      digiCfg.outputMeasurementParticlesMap = "measurement_particles_map";
      digiCfg.outputMeasurementSimHitsMap = "measurement_simhits_map";
      digiCfg.trackingGeometry = trackingGeometry;
      digiCfg.randomNumbers = randomNumbers;
      digiCfg.digitizationConfigs = digitizationConfigs;

      sequencer.addAlgorithm(
          std::make_shared<DigitizationAlgorithm>(digiCfg, logLevel));

      // Write digitization output as ROOT files
      if (vars["output-root"].template as<bool>()) {
        RootDigitizationWriter::Config digitWriterRoot;
        digitWriterRoot.inputMeasurements = digiCfg.outputMeasurements;
        digitWriterRoot.inputSimHits = digiCfg.inputSimHits;
        digitWriterRoot.inputMeasurementSimHitsMap =
            digiCfg.outputMeasurementSimHitsMap;
        digitWriterRoot.filePath =
            joinPaths(outputDir, digiCfg.outputMeasurements + ".root");
        // Translate into bound indices of the digitization configuration
        std::vector<std::pair<Acts::GeometryIdentifier,
                              std::vector<Acts::BoundIndices>>>
            bIndexInput;
        for (size_t ibi = 0; ibi < digiCfg.digitizationConfigs.size(); ++ibi) {
          Acts::GeometryIdentifier geoID =
              digiCfg.digitizationConfigs.idAt(ibi);
          const auto dCfg = digiCfg.digitizationConfigs.valueAt(ibi);
          std::vector<Acts::BoundIndices> boundIndices;
          boundIndices.insert(boundIndices.end(),
                              dCfg.geometricDigiConfig.indices.begin(),
                              dCfg.geometricDigiConfig.indices.end());
          for (const auto& sConfig : dCfg.smearingDigiConfig) {
            boundIndices.push_back(sConfig.index);
          }
          bIndexInput.push_back({geoID, boundIndices});
        }
        digitWriterRoot.boundIndices =
            Acts::GeometryHierarchyMap<std::vector<Acts::BoundIndices>>(
                bIndexInput);

        digitWriterRoot.trackingGeometry = trackingGeometry;
        sequencer.addWriter(std::make_shared<RootDigitizationWriter>(
            digitWriterRoot, logLevel));
      }
    }

  } else if (vars["digi-smear"].as<bool>()) {
    SmearingAlgorithm::Config smearCfg = Options::readSmearingConfig(vars);
    smearCfg.inputSimHits = kFatrasCollectionHits;
    smearCfg.outputMeasurements = "measurements";
    smearCfg.outputSourceLinks = "sourcelinks";
    smearCfg.outputMeasurementParticlesMap = "measurement_particles_map";
    smearCfg.outputMeasurementSimHitsMap = "measurement_simhits_map";
    smearCfg.trackingGeometry = trackingGeometry;
    smearCfg.randomNumbers = randomNumbers;
    sequencer.addAlgorithm(
        std::make_shared<SmearingAlgorithm>(smearCfg, logLevel));

    // Write digitization output as ROOT files
    if (vars["output-root"].template as<bool>()) {
      RootDigitizationWriter::Config smearWriterRoot;
      smearWriterRoot.inputMeasurements = smearCfg.outputMeasurements;
      smearWriterRoot.inputSimHits = smearCfg.inputSimHits;
      smearWriterRoot.inputMeasurementSimHitsMap =
          smearCfg.outputMeasurementSimHitsMap;
      smearWriterRoot.filePath =
          joinPaths(outputDir, smearCfg.outputMeasurements + ".root");
      // Translate into bound indices of the digitization configuration
      std::vector<
          std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
          bIndexInput;
      for (size_t ibi = 0; ibi < smearCfg.smearers.size(); ++ibi) {
        Acts::GeometryIdentifier geoID = smearCfg.smearers.idAt(ibi);
        const auto sCfg = smearCfg.smearers.valueAt(ibi);
        std::vector<Acts::BoundIndices> boundIndices;

        for (const auto& sConfig : sCfg) {
          boundIndices.push_back(sConfig.index);
        }
        bIndexInput.push_back({geoID, boundIndices});
      }
      smearWriterRoot.boundIndices =
          Acts::GeometryHierarchyMap<std::vector<Acts::BoundIndices>>(
              bIndexInput);
      smearWriterRoot.trackingGeometry = trackingGeometry;
      sequencer.addWriter(
          std::make_shared<RootDigitizationWriter>(smearWriterRoot, logLevel));
    }

  } else if (vars["digi-geometric-3d"].as<bool>()) {
    // Configure the digitizer
    PlanarSteppingAlgorithm::Config digi;
    digi.inputSimHits = "hits";
    digi.outputClusters = "clusters";
    digi.outputMeasurements = "measurements";
    digi.outputSourceLinks = "sourcelinks";
    digi.outputMeasurementParticlesMap = "measurement_particles_map";
    digi.outputMeasurementSimHitsMap = "measurement_simhits_map";
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
      clusterWriterCsv.inputSimHits = digi.inputSimHits;
      clusterWriterCsv.outputDir = outputDir;
      clusterWriterCsv.trackingGeometry = trackingGeometry;
      sequencer.addWriter(
          std::make_shared<CsvPlanarClusterWriter>(clusterWriterCsv, logLevel));
    }

    // Write digitization output as ROOT files
    if (vars["output-root"].template as<bool>()) {
      // clusters as root
      RootPlanarClusterWriter::Config clusterWriterRoot;
      clusterWriterRoot.inputClusters = digi.outputClusters;
      clusterWriterRoot.inputSimHits = digi.inputSimHits;
      clusterWriterRoot.filePath =
          joinPaths(outputDir, digi.outputClusters + ".root");
      sequencer.addWriter(std::make_shared<RootPlanarClusterWriter>(
          clusterWriterRoot, logLevel));
    }
  }
}
