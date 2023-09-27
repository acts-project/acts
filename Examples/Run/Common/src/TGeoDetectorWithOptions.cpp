// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/TGeoDetectorWithOptions.hpp"

#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Plugins/TGeo/TGeoCylinderDiscSplitter.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/TGeoDetector/JsonTGeoDetectorConfig.hpp"
#include "ActsExamples/TGeoDetector/TGeoITkModuleSplitter.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <cstdlib>
#include <fstream>
#include <list>

#include <boost/program_options.hpp>

namespace ActsExamples {
using namespace Options;

namespace {

/// Read the TGeo layer builder configurations from the user configuration
/// specified with --geo-tgeo-jsonconfig.
void readTGeoLayerBuilderConfigs(const Variables& vm,
                                 TGeoDetector::Config& config) {
  const auto path = vm["geo-tgeo-jsonconfig"].template as<std::string>();
  TGeoDetector::readTGeoLayerBuilderConfigsFile(path, config);
}

/// Dump TGeo Detector config to file.
void writeTGeoDetectorConfig(const Variables& vm,
                             TGeoDetector::Config& config) {
  const auto path = vm["geo-tgeo-dump-jsonconfig"].template as<std::string>();
  nlohmann::json djson;
  if (path.empty()) {
    return;
  }
  std::ofstream outfile(path, std::ofstream::out | std::ofstream::binary);

  djson["geo-tgeo-unit-scalor"] = config.unitScalor;
  djson["geo-tgeo-build-beampipe"] = config.buildBeamPipe;
  djson["geo-tgeo-beampipe-parameters"] =
      std::array<double, 3>{config.beamPipeRadius, config.beamPipeHalflengthZ,
                            config.beamPipeLayerThickness};

  // Enable empty volume dump
  if (config.volumes.empty()) {
    config.volumes.emplace_back();
  }
  djson["Volumes"] = config.volumes;

  outfile << djson.dump(2) << std::endl;
}

}  // namespace

void TGeoDetectorWithOptions::addOptions(
    boost::program_options::options_description& opt) const {
  using boost::program_options::value;

  auto tmp = opt.add_options();
  // required global options
  tmp("geo-tgeo-filename", value<std::string>()->default_value(""),
      "Root file name.");
  tmp("geo-tgeo-jsonconfig", value<std::string>()->default_value(""),
      "Json config file name.");
  tmp("geo-tgeo-dump-jsonconfig",
      value<std::string>()->default_value("tgeo_empty_config.json"),
      "Json file to dump empty config into.");
}

auto TGeoDetectorWithOptions::finalize(
    const boost::program_options::variables_map& vm,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  TGeoDetector::Config config;

  config.fileName = vm["geo-tgeo-filename"].as<std::string>();

  config.surfaceLogLevel =
      Acts::Logging::Level(vm["geo-surface-loglevel"].template as<size_t>());
  config.layerLogLevel =
      Acts::Logging::Level(vm["geo-layer-loglevel"].template as<size_t>());
  config.volumeLogLevel =
      Acts::Logging::Level(vm["geo-volume-loglevel"].template as<size_t>());

  // No valid geometry configuration. Stop
  if (vm["geo-tgeo-jsonconfig"].as<std::string>().empty()) {
    writeTGeoDetectorConfig(vm, config);
    std::exit(EXIT_SUCCESS);
  }
  // Enable dump from full config
  else if (not(vm["geo-tgeo-dump-jsonconfig"].as<std::string>().compare(
                   "tgeo_empty_cofig.json") == 0)) {
    readTGeoLayerBuilderConfigs(vm, config);
    writeTGeoDetectorConfig(vm, config);
  } else {
    readTGeoLayerBuilderConfigs(vm, config);
  }

  return m_detector.finalize(config, std::move(mdecorator));
}

}  // namespace ActsExamples
