// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "ActsExamples/Digitization/SmearingAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsWriter.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <fstream>
#include <memory>

#include <boost/program_options.hpp>

#include "DigitizationInput.hpp"

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int runDigitizationExample(
    int argc, char* argv[],
    std::shared_ptr<ActsExamples::IBaseDetector> detector) {
  // Setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::Csv | OutputFormat::Root);
  Options::addCsvWriterOptions(desc);
  Options::addInputOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addDigitizationOptions(desc);

  auto opt = desc.add_options();
  opt("digi-read-write-test", boost::program_options::bool_switch(),
      "Test reading and writing roundtrip.");

  // Add specific options for this geometry
  detector->addOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = Options::readLogLevel(vm);
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());

  // The geometry, material and decoration
  auto geometry = Geometry::build(vm, *detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  auto randomNumbers =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vm));

  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Setup the magnetic field
  Options::setupMagneticFieldServices(vm, sequencer);
  auto magneticField = Options::readMagneticField(vm);

  // Read the sim hits
  auto simHitReaderCfg = setupSimHitReading(vm, sequencer);
  // Read the particles
  auto particleReaderCfg = setupParticleReading(vm, sequencer);

  auto cfile = vm["digi-config-file"].as<std::string>();

  if (not cfile.empty() or vm["digi-smear"].as<bool>()) {
    // To be handled differently between DigitizationAlgorithm/SmearingAlgorithm
    std::string outputMeasurements = "measurements";
    std::string outputSourceLinks = "sourcelinks";
    std::string outputClusters = "";
    std::string outputMeasurementParticlesMap = "measurement_particles_map";
    std::string outputMeasurementSimHitsMap = "measurement_simhits_map";

    std::vector<
        std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
        bIndexInput;

    if (not cfile.empty()) {
      // Configuration by json file evokes 'DigitizationAlgorithm'
      auto in = std::ifstream(cfile, std::ifstream::in | std::ifstream::binary);
      if (in.good()) {
        // Get the json file for the configuration
        nlohmann::json djson;
        in >> djson;
        in.close();

        outputClusters = "clusters";

        DigiConfigContainer digitizationConfigs =
            DigiConfigConverter("digitization-configuration").fromJson(djson);

        if (vm["digi-read-write-test"].as<bool>()) {
          // Read measurements from CSV files
          auto measReaderCfg = Options::readCsvMeasurementReaderConfig(vm);
          measReaderCfg.inputDir = vm["input-dir"].as<std::string>();
          measReaderCfg.outputMeasurements = outputMeasurements;
          measReaderCfg.outputClusters = outputClusters;
          sequencer.addReader(
              std::make_shared<CsvMeasurementReader>(measReaderCfg, logLevel));

        } else {
          DigitizationAlgorithm::Config digiCfg =
              Options::readDigitizationConfig(vm);
          digiCfg.inputSimHits = simHitReaderCfg.outputSimHits;
          digiCfg.outputMeasurements = outputMeasurements;
          digiCfg.outputClusters = outputClusters;
          digiCfg.outputSourceLinks = outputSourceLinks;
          digiCfg.outputMeasurementParticlesMap = outputMeasurementSimHitsMap;
          digiCfg.outputMeasurementSimHitsMap = outputMeasurementParticlesMap;
          digiCfg.trackingGeometry = tGeometry;
          digiCfg.randomNumbers = randomNumbers;
          digiCfg.digitizationConfigs = digitizationConfigs;
          sequencer.addAlgorithm(
              std::make_shared<DigitizationAlgorithm>(digiCfg, logLevel));
        }

        // Prepare for the (eventual) output writing
        for (size_t ibi = 0; ibi < digitizationConfigs.size(); ++ibi) {
          Acts::GeometryIdentifier geoID = digitizationConfigs.idAt(ibi);
          const auto dCfg = digitizationConfigs.valueAt(ibi);
          std::vector<Acts::BoundIndices> boundIndices;
          boundIndices.insert(boundIndices.end(),
                              dCfg.geometricDigiConfig.indices.begin(),
                              dCfg.geometricDigiConfig.indices.end());
          for (const auto& sConfig : dCfg.smearingDigiConfig) {
            boundIndices.push_back(sConfig.index);
          }
          bIndexInput.push_back({geoID, boundIndices});
        }
      }
    } else if (vm["digi-smear"].as<bool>()) {
      // Simpler configuration for smearing algorithm
      SmearingAlgorithm::Config smearCfg = Options::readSmearingConfig(vm);
      smearCfg.inputSimHits = simHitReaderCfg.outputSimHits;
      smearCfg.outputMeasurements = outputMeasurements;
      smearCfg.outputSourceLinks = outputSourceLinks;
      smearCfg.outputMeasurementParticlesMap = outputMeasurementParticlesMap;
      smearCfg.outputMeasurementSimHitsMap = outputMeasurementSimHitsMap;
      smearCfg.trackingGeometry = tGeometry;
      smearCfg.randomNumbers = randomNumbers;
      sequencer.addAlgorithm(
          std::make_shared<SmearingAlgorithm>(smearCfg, logLevel));

      // Prepare for the (eventual) output writing
      for (size_t ibi = 0; ibi < smearCfg.smearers.size(); ++ibi) {
        Acts::GeometryIdentifier geoID = smearCfg.smearers.idAt(ibi);
        const auto sCfg = smearCfg.smearers.valueAt(ibi);
        std::vector<Acts::BoundIndices> boundIndices;

        for (const auto& sConfig : sCfg) {
          boundIndices.push_back(sConfig.index);
        }
        bIndexInput.push_back({geoID, boundIndices});
      }
    }

    // Write digitization output as ROOT files
    if (vm["output-root"].template as<bool>()) {
      RootMeasurementWriter::Config measWriterRoot;
      measWriterRoot.inputMeasurements = outputMeasurements;
      measWriterRoot.inputClusters = outputClusters;
      measWriterRoot.inputSimHits = simHitReaderCfg.outputSimHits;
      measWriterRoot.inputMeasurementSimHitsMap = outputMeasurementSimHitsMap;
      measWriterRoot.filePath =
          joinPaths(outputDir, std::string(outputMeasurements) + ".root");
      measWriterRoot.boundIndices =
          Acts::GeometryHierarchyMap<std::vector<Acts::BoundIndices>>(
              bIndexInput);
      measWriterRoot.trackingGeometry = tGeometry;
      sequencer.addWriter(
          std::make_shared<RootMeasurementWriter>(measWriterRoot, logLevel));
    }

    // Write digitization out as CSV files
    if (vm["output-csv"].template as<bool>()) {
      CsvMeasurementWriter::Config measWriterCsv =
          Options::readCsvMeasurementWriterConfig(vm);
      measWriterCsv.inputMeasurements = outputMeasurements;
      measWriterCsv.inputClusters = outputClusters;
      measWriterCsv.inputSimHits = simHitReaderCfg.outputSimHits;

      measWriterCsv.inputMeasurementSimHitsMap = outputMeasurementSimHitsMap;
      sequencer.addWriter(
          std::make_shared<CsvMeasurementWriter>(measWriterCsv, logLevel));
    }
  }

  return sequencer.run();
}
