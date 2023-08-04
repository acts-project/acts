// This file is part of the Acts project.
//
// Copyright (C) 2018-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvTrackingGeometryWriter.hpp"
#include "ActsExamples/Io/Json/JsonMaterialWriter.hpp"
#include "ActsExamples/Io/Json/JsonSurfacesWriter.hpp"
#include "ActsExamples/Io/Root/RootMaterialWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/CsvOptionsWriter.hpp"
#include "ActsExamples/Options/JsonOptionsWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjTrackingGeometryWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjWriterOptions.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <string>
#include <vector>

int processGeometry(int argc, char* argv[],
                    ActsExamples::IBaseDetector& detector) {
  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addGeometryOptions(desc);
  ActsExamples::Options::addMaterialOptions(desc);
  ActsExamples::Options::addObjWriterOptions(desc);
  ActsExamples::Options::addCsvWriterOptions(desc);
  ActsExamples::Options::addJsonWriterOptions(desc);
  ActsExamples::Options::addOutputOptions(
      desc,
      ActsExamples::OutputFormat::Root | ActsExamples::OutputFormat::Json |
          ActsExamples::OutputFormat::Cbor | ActsExamples::OutputFormat::Csv |
          ActsExamples::OutputFormat::Obj);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  // Now read the standard options
  auto logLevel = ActsExamples::Options::readLogLevel(vm);
  size_t nEvents =
      ActsExamples::Options::readSequencerConfig(vm).events.value_or(1);

  // The geometry, material and decoration
  auto geometry = ActsExamples::Geometry::build(vm, detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;

  // The detectors
  auto volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].as<size_t>());

  for (size_t ievt = 0; ievt < nEvents; ++ievt) {
    // Setup the event and algorithm context
    ActsExamples::WhiteBoard eventStore(
        Acts::getDefaultLogger("EventStore#" + std::to_string(ievt), logLevel));
    size_t ialg = 0;

    // The geometry context
    ActsExamples::AlgorithmContext context(ialg, ievt, eventStore);

    /// Decorate the context
    for (auto& cdr : contextDecorators) {
      if (cdr->decorate(context) != ActsExamples::ProcessCode::SUCCESS) {
        throw std::runtime_error("Failed to decorate event context");
      }
    }

    std::string geoContextStr = "";
    if (!contextDecorators.empty()) {
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
      auto tgObjWriterConfig =
          ActsExamples::Options::readObjTrackingGeometryWriterConfig(vm);
      tgObjWriterConfig.outputDir = outputDir;
      auto tgObjWriter =
          std::make_shared<ActsExamples::ObjTrackingGeometryWriter>(
              tgObjWriterConfig, volumeLogLevel);
      // Write the tracking geometry object
      tgObjWriter->write(context, *tGeometry);
    }

    // CSV output
    if (vm["output-csv"].as<bool>()) {
      // setup the tracking geometry writer
      ActsExamples::CsvTrackingGeometryWriter::Config tgCsvWriterConfig;
      tgCsvWriterConfig.trackingGeometry = tGeometry;
      tgCsvWriterConfig.outputDir = outputDir;
      tgCsvWriterConfig.writePerEvent = true;
      auto tgCsvWriter =
          std::make_shared<ActsExamples::CsvTrackingGeometryWriter>(
              tgCsvWriterConfig, logLevel);

      // Write the tracking geometry object
      tgCsvWriter->write(context);
    }

    // JSON output
    if (vm["output-json"].as<bool>()) {
      auto sJsonWriterConfig =
          ActsExamples::Options::readJsonSurfacesWriterConfig(vm);
      sJsonWriterConfig.trackingGeometry = tGeometry;
      sJsonWriterConfig.outputDir = outputDir;
      sJsonWriterConfig.writePerEvent = true;
      auto sJsonWriter = std::make_shared<ActsExamples::JsonSurfacesWriter>(
          sJsonWriterConfig, logLevel);

      // Write the tracking geometry object
      sJsonWriter->write(context);
    }

    // Get the file name from the options
    std::string materialFileName = vm["mat-output-file"].as<std::string>();

    if (!materialFileName.empty() and vm["output-root"].template as<bool>()) {
      // The writer of the indexed material
      ActsExamples::RootMaterialWriter::Config rmwConfig;
      rmwConfig.filePath = materialFileName + ".root";
      ActsExamples::RootMaterialWriter rmwImpl(rmwConfig, logLevel);
      rmwImpl.write(*tGeometry);
    }

    if (!materialFileName.empty() and (vm["output-json"].template as<bool>() or
                                       vm["output-cbor"].template as<bool>())) {
      // the material writer
      Acts::MaterialMapJsonConverter::Config jmConverterCfg;
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
      jmConverterCfg.processNonMaterial =
          vm["mat-output-allmaterial"].template as<bool>();
      jmConverterCfg.context = context.geoContext;
      // The writer
      ActsExamples::JsonMaterialWriter::Config jmWriterCfg;
      jmWriterCfg.converterCfg = std::move(jmConverterCfg);
      jmWriterCfg.fileName = materialFileName;
      ActsExamples::JsonFormat format = ActsExamples::JsonFormat::NoOutput;
      if (vm["output-json"].template as<bool>()) {
        format = format | ActsExamples::JsonFormat::Json;
      }
      if (vm["output-cbor"].template as<bool>()) {
        format = format | ActsExamples::JsonFormat::Cbor;
      }
      jmWriterCfg.writeFormat = format;

      ActsExamples::JsonMaterialWriter jmwImpl(jmWriterCfg, logLevel);

      jmwImpl.write(*tGeometry);
    }
  }

  return 0;
}
