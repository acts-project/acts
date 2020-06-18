// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <memory>
#include <string>
#include <vector>

#include "ACTFW/Detector/IBaseDetector.hpp"
#include "ACTFW/Framework/AlgorithmContext.hpp"
#include "ACTFW/Framework/IContextDecorator.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Csv/CsvOptionsWriter.hpp"
#include "ACTFW/Io/Csv/CsvTrackingGeometryWriter.hpp"
#include "ACTFW/Io/Root/RootMaterialWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/Json/JsonMaterialWriter.hpp"
#include "ACTFW/Plugins/Obj/ObjTrackingGeometryWriter.hpp"
#include "ACTFW/Plugins/Obj/ObjWriterOptions.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"

int processGeometry(int argc, char* argv[], FW::IBaseDetector& detector) {
  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addObjWriterOptions(desc);
  FW::Options::addCsvWriterOptions(desc);
  FW::Options::addOutputOptions(desc);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // Now read the standard options
  auto logLevel = FW::Options::readLogLevel(vm);
  auto nEvents = FW::Options::readSequencerConfig(vm).events;

  // The geometry, material and decoration
  auto geometry = FW::Geometry::build(vm, detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;

  // The detectors
  read_strings subDetectors = vm["geo-detector-volume"].as<read_strings>();

  auto volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].as<size_t>());

  for (size_t ievt = 0; ievt < nEvents; ++ievt) {
    // Setup the event and algorithm context
    FW::WhiteBoard eventStore(
        Acts::getDefaultLogger("EventStore#" + std::to_string(ievt), logLevel));
    size_t ialg = 0;

    // The geometry context
    FW::AlgorithmContext context(ialg, ievt, eventStore);

    /// Decorate the context
    for (auto& cdr : contextDecorators) {
      if (cdr->decorate(context) != FW::ProcessCode::SUCCESS)
        throw std::runtime_error("Failed to decorate event context");
    }

    std::string geoContextStr = "";
    if (contextDecorators.size() > 0) {
      // We need indeed a context object
      if (nEvents > 1) {
        geoContextStr = "_geoContext" + std::to_string(ievt);
      }
    }

    // ---------------------------------------------------------------------------------
    // Output directory
    std::string outputDir = vm["output-dir"].template as<std::string>();

    // OBJ output
    if (vm["output-obj"].as<bool>()) {
      // Configure the tracking geometry writer
      auto tgObjWriterConfig = FW::Options::readObjTrackingGeometryWriterConfig(
          vm, "ObjTrackingGeometryWriter", volumeLogLevel);
      auto tgObjWriter = std::make_shared<FW::Obj::ObjTrackingGeometryWriter>(
          tgObjWriterConfig);
      // Write the tracking geometry object
      tgObjWriter->write(context, *tGeometry);
    }

    // CSV output
    if (vm["output-csv"].as<bool>()) {
      // setup the tracking geometry writer
      FW::CsvTrackingGeometryWriter::Config tgCsvWriterConfig;
      tgCsvWriterConfig.trackingGeometry = tGeometry;
      tgCsvWriterConfig.outputDir = outputDir;
      tgCsvWriterConfig.writePerEvent = true;
      auto tgCsvWriter = std::make_shared<FW::CsvTrackingGeometryWriter>(
          tgCsvWriterConfig, logLevel);

      // Write the tracking geometry object
      tgCsvWriter->write(context);
    }

    // Get the file name from the options
    std::string materialFileName = vm["mat-output-file"].as<std::string>();

    if (!materialFileName.empty() and vm["output-root"].template as<bool>()) {
      // The writer of the indexed material
      FW::RootMaterialWriter::Config rmwConfig("MaterialWriter");
      rmwConfig.fileName = materialFileName + ".root";
      FW::RootMaterialWriter rmwImpl(rmwConfig);
      rmwImpl.write(*tGeometry);
    }

    if (!materialFileName.empty() and vm["output-json"].template as<bool>()) {
      /// The name of the output file
      std::string fileName = vm["mat-output-file"].template as<std::string>();
      // the material writer
      Acts::JsonGeometryConverter::Config jmConverterCfg(
          "JsonGeometryConverter", Acts::Logging::INFO);
      jmConverterCfg.processSensitives =
          vm["mat-output-sensitives"].template as<bool>();
      jmConverterCfg.processApproaches =
          vm["mat-output-approaches"].template as<bool>();
      jmConverterCfg.processRepresenting =
          vm["mat-output-representing"].template as<bool>();
      jmConverterCfg.processBoundaries =
          vm["mat-output-boundaries"].template as<bool>();
      jmConverterCfg.processVolumes =
          vm["mat-output-volumes"].template as<bool>();
      jmConverterCfg.processDenseVolumes =
          vm["mat-output-dense-volumes"].template as<bool>();
      jmConverterCfg.writeData = vm["mat-output-data"].template as<bool>();
      jmConverterCfg.processnonmaterial =
          vm["mat-output-allmaterial"].template as<bool>();
      // The writer
      FW::Json::JsonMaterialWriter jmwImpl(std::move(jmConverterCfg),
                                           materialFileName + ".json");

      jmwImpl.write(*tGeometry);
    }
  }

  return 0;
}
