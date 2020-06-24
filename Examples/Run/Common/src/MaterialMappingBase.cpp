// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Material/SurfaceMaterialMapper.hpp>
#include <Acts/Plugins/Json/JsonGeometryConverter.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/StraightLineStepper.hpp>
#include <boost/program_options.hpp>
#include <memory>

#include "ACTFW/Detector/IBaseDetector.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Root/RootMaterialTrackReader.hpp"
#include "ACTFW/Io/Root/RootMaterialTrackWriter.hpp"
#include "ACTFW/Io/Root/RootMaterialWriter.hpp"
#include "ACTFW/MaterialMapping/MaterialMapping.hpp"
#include "ACTFW/MaterialMapping/MaterialMappingOptions.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/Json/JsonMaterialWriter.hpp"
#include "ACTFW/Propagation/PropagationOptions.hpp"
#include "ACTFW/Utilities/Paths.hpp"

namespace po = boost::program_options;

int materialMappingExample(int argc, char* argv[],
                           FW::IBaseDetector& detector) {
  // Setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addMaterialMappingOptions(desc);
  FW::Options::addPropagationOptions(desc);
  FW::Options::addInputOptions(desc);
  FW::Options::addOutputOptions(desc);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  FW::Sequencer sequencer(FW::Options::readSequencerConfig(vm));

  // Get the log level
  auto logLevel = FW::Options::readLogLevel(vm);

  // The geometry, material and decoration
  auto geometry = FW::Geometry::build(vm, detector);
  auto tGeometry = geometry.first;

  /// Default contexts
  Acts::GeometryContext geoContext;
  Acts::MagneticFieldContext mfContext;

  // Straight line stepper
  using SlStepper = Acts::StraightLineStepper;
  using Propagator = Acts::Propagator<SlStepper, Acts::Navigator>;

  auto matCollection = vm["mat-mapping-collection"].template as<std::string>();
  auto mapSurface = vm["mat-mapping-surfaces"].template as<bool>();
  auto mapVolume = vm["mat-mapping-volumes"].template as<bool>();
  auto volumeStep = vm["mat-mapping-volume-stepsize"].template as<float>();
  if (!mapSurface && !mapVolume) {
    return EXIT_FAILURE;
  }
  // ---------------------------------------------------------------------------------
  // Input directory & input file handling
  std::string intputDir = vm["input-dir"].template as<std::string>();
  auto intputFiles = vm["input-files"].template as<read_strings>();

  if (vm["input-root"].template as<bool>()) {
    // Read the material step information from a ROOT TTree
    FW::RootMaterialTrackReader::Config matTrackReaderRootConfig;
    if (not matCollection.empty()) {
      matTrackReaderRootConfig.collection = matCollection;
    }
    matTrackReaderRootConfig.fileList = intputFiles;
    auto matTrackReaderRoot =
        std::make_shared<FW::RootMaterialTrackReader>(matTrackReaderRootConfig);
    sequencer.addReader(matTrackReaderRoot);
  }

  /// The material mapping algorithm
  FW::MaterialMapping::Config mmAlgConfig(geoContext, mfContext);
  if (mapSurface) {
    // Get a Navigator
    Acts::Navigator navigator(tGeometry);
    // Make stepper and propagator
    SlStepper stepper;
    Propagator propagator(std::move(stepper), std::move(navigator));
    /// The material surface mapper
    Acts::SurfaceMaterialMapper::Config smmConfig;
    auto smm = std::make_shared<Acts::SurfaceMaterialMapper>(
        smmConfig, std::move(propagator),
        Acts::getDefaultLogger("SurfaceMaterialMapper", logLevel));
    mmAlgConfig.materialSurfaceMapper = smm;
  }
  if (mapVolume) {
    // Get a Navigator
    Acts::Navigator navigator(tGeometry);
    // Make stepper and propagator
    SlStepper stepper;
    Propagator propagator(std::move(stepper), std::move(navigator));
    /// The material volume mapper
    Acts::VolumeMaterialMapper::Config vmmConfig;
    vmmConfig.mappingStep = volumeStep;
    auto vmm = std::make_shared<Acts::VolumeMaterialMapper>(
        vmmConfig, std::move(propagator),
        Acts::getDefaultLogger("VolumeMaterialMapper", logLevel));
    mmAlgConfig.materialVolumeMapper = vmm;
  }
  mmAlgConfig.trackingGeometry = tGeometry;

  // Get the file name from the options
  std::string materialFileName = vm["mat-output-file"].as<std::string>();

  if (!materialFileName.empty() and vm["output-root"].template as<bool>()) {
    // The writer of the indexed material
    FW::RootMaterialWriter::Config rmwConfig("MaterialWriter");
    rmwConfig.fileName = materialFileName + ".root";
    FW::RootMaterialWriter rmwImpl(rmwConfig);
    // Fullfill the IMaterialWriter interface
    using RootWriter = FW::MaterialWriterT<FW::RootMaterialWriter>;
    mmAlgConfig.materialWriters.push_back(
        std::make_shared<RootWriter>(std::move(rmwImpl)));

    if (mapSurface) {
      // Write the propagation steps as ROOT TTree
      FW::RootMaterialTrackWriter::Config matTrackWriterRootConfig;
      matTrackWriterRootConfig.filePath = materialFileName + "_tracks.root";
      matTrackWriterRootConfig.collection =
          mmAlgConfig.mappingMaterialCollection;
      matTrackWriterRootConfig.storesurface = true;
      auto matTrackWriterRoot = std::make_shared<FW::RootMaterialTrackWriter>(
          matTrackWriterRootConfig, logLevel);
      sequencer.addWriter(matTrackWriterRoot);
    }
  }

  if (!materialFileName.empty() and vm["output-json"].template as<bool>()) {
    /// The name of the output file
    std::string fileName = vm["mat-output-file"].template as<std::string>();
    // the material writer
    Acts::JsonGeometryConverter::Config jmConverterCfg("JsonGeometryConverter",
                                                       Acts::Logging::INFO);
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
    jmConverterCfg.writeData = vm["mat-output-data"].template as<bool>();
    // The writer
    FW::Json::JsonMaterialWriter jmwImpl(jmConverterCfg,
                                         materialFileName + ".json");
    // Fullfill the IMaterialWriter interface
    using JsonWriter = FW::MaterialWriterT<FW::Json::JsonMaterialWriter>;
    mmAlgConfig.materialWriters.push_back(
        std::make_shared<JsonWriter>(std::move(jmwImpl)));
  }

  // Create the material mapping
  auto mmAlg = std::make_shared<FW::MaterialMapping>(mmAlgConfig);

  // Append the Algorithm
  sequencer.addAlgorithm(mmAlg);

  // Initiate the run
  sequencer.run();
  // Return success code
  return 0;
}
