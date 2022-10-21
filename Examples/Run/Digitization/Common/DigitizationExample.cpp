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
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/CsvOptionsReader.hpp"
#include "ActsExamples/Options/CsvOptionsWriter.hpp"
#include "ActsExamples/Options/DigitizationOptions.hpp"
#include "ActsExamples/Options/MagneticFieldOptions.hpp"
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

  auto digiCfg = DigitizationConfig(
      vm["digi-merge"].as<bool>(), vm["digi-merge-nsigma"].as<double>(),
      vm["digi-merge-common-corner"].as<bool>(),
      readDigiConfigFromJson(vm["digi-config-file"].as<std::string>()));
  digiCfg.inputSimHits = simHitReaderCfg.outputSimHits;
  digiCfg.trackingGeometry = tGeometry;
  digiCfg.randomNumbers = randomNumbers;

  if (not vm["dump-digi-config"].as<std::string>().empty()) {
    writeDigiConfigToJson(digiCfg.digitizationConfigs,
                          vm["dump-digi-config"].as<std::string>());
  }

  std::vector<
      std::pair<Acts::GeometryIdentifier, std::vector<Acts::BoundIndices>>>
      bIndexInput;

  if (vm["digi-read-write-test"].as<bool>()) {
    // Read measurements from CSV files
    auto measReaderCfg = Options::readCsvMeasurementReaderConfig(vm);
    measReaderCfg.inputDir = vm["input-dir"].as<std::string>();
    measReaderCfg.outputMeasurements = digiCfg.outputMeasurements;
    measReaderCfg.outputSourceLinks = digiCfg.outputSourceLinks;
    measReaderCfg.outputClusters = digiCfg.outputClusters;
    measReaderCfg.outputMeasurementSimHitsMap =
        digiCfg.outputMeasurementSimHitsMap;
    sequencer.addReader(
        std::make_shared<CsvMeasurementReader>(measReaderCfg, logLevel));
  } else {
    sequencer.addAlgorithm(
        std::make_shared<DigitizationAlgorithm>(digiCfg, logLevel));
  }

  // Write digitization output as ROOT files
  if (vm["output-root"].template as<bool>()) {
    RootMeasurementWriter::Config measWriterRoot;
    measWriterRoot.inputMeasurements = digiCfg.outputMeasurements;
    measWriterRoot.inputClusters = digiCfg.outputClusters;
    measWriterRoot.inputSimHits = simHitReaderCfg.outputSimHits;
    measWriterRoot.inputMeasurementSimHitsMap =
        digiCfg.outputMeasurementSimHitsMap;
    measWriterRoot.filePath =
        joinPaths(outputDir, std::string(digiCfg.outputMeasurements) + ".root");
    measWriterRoot.boundIndices =
        Acts::GeometryHierarchyMap<std::vector<Acts::BoundIndices>>(
            digiCfg.getBoundIndices());
    measWriterRoot.trackingGeometry = tGeometry;
    sequencer.addWriter(
        std::make_shared<RootMeasurementWriter>(measWriterRoot, logLevel));
  }

  // Write digitization out as CSV files
  if (vm["output-csv"].template as<bool>()) {
    CsvMeasurementWriter::Config measWriterCsv =
        Options::readCsvMeasurementWriterConfig(vm);
    measWriterCsv.inputMeasurements = digiCfg.outputMeasurements;
    measWriterCsv.inputClusters = digiCfg.outputClusters;
    measWriterCsv.inputSimHits = simHitReaderCfg.outputSimHits;
    measWriterCsv.inputMeasurementSimHitsMap =
        digiCfg.outputMeasurementSimHitsMap;
    sequencer.addWriter(
        std::make_shared<CsvMeasurementWriter>(measWriterCsv, logLevel));
  }

  return sequencer.run();
}
