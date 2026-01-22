// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"

#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ITrackingVolumeBuilder.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/ProtoLayerHelper.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/SurfaceBinningMatcher.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/TGeoDetector/JsonTGeoDetectorConfig.hpp"
#include "ActsExamples/TGeoDetector/TGeoITkModuleSplitter.hpp"
#include "ActsPlugins/Root/TGeoCylinderDiscSplitter.hpp"
#include "ActsPlugins/Root/TGeoLayerBuilder.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <limits>
#include <list>
#include <optional>
#include <utility>

#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>

#include "TGeoManager.h"

using namespace Acts;
using namespace ActsPlugins;

namespace ActsExamples {

namespace {

/// @brief Function that constructs a set of layer builder
///        configs from a central @c TGeoDetector config.
///
/// @param config The input config
/// @return Vector of layer builder configs
std::vector<ActsPlugins::TGeoLayerBuilder::Config> makeLayerBuilderConfigs(
    const TGeoDetector::Config& config, const Logger& logger) {
  std::vector<ActsPlugins::TGeoLayerBuilder::Config> detLayerConfigs;

  // iterate over all configured detector volumes
  for (const auto& volume : config.volumes) {
    ActsPlugins::TGeoLayerBuilder::Config layerBuilderConfig;
    layerBuilderConfig.configurationName = volume.name;
    layerBuilderConfig.unit = config.unitScalor;
    layerBuilderConfig.detectorElementFactory = config.detectorElementFactory;

    // configure surface autobinning
    std::vector<std::pair<double, double>> binTolerances(numAxisDirections(),
                                                         {0., 0.});
    binTolerances[toUnderlying(AxisDirection::AxisR)] = {
        volume.binToleranceR.lower.value_or(0.),
        volume.binToleranceR.upper.value_or(0.)};
    binTolerances[toUnderlying(AxisDirection::AxisZ)] = {
        volume.binToleranceZ.lower.value_or(0.),
        volume.binToleranceZ.upper.value_or(0.)};
    binTolerances[toUnderlying(AxisDirection::AxisPhi)] = {
        volume.binTolerancePhi.lower.value_or(0.),
        volume.binTolerancePhi.upper.value_or(0.)};

    layerBuilderConfig.autoSurfaceBinning = true;
    layerBuilderConfig.surfaceBinMatcher = SurfaceBinningMatcher(binTolerances);

    // loop over the negative/central/positive layer configurations
    for (auto ncp : {
             TGeoDetector::Config::Negative,
             TGeoDetector::Config::Central,
             TGeoDetector::Config::Positive,
         }) {
      if (!volume.layers.at(ncp)) {
        continue;
      }

      ActsPlugins::TGeoLayerBuilder::LayerConfig lConfig;
      lConfig.volumeName = volume.subVolumeName.at(ncp);
      lConfig.sensorNames = volume.sensitiveNames.at(ncp);
      lConfig.localAxes = volume.sensitiveAxes.at(ncp);
      lConfig.envelope = {config.layerEnvelopeR, config.layerEnvelopeR};

      auto rR = volume.rRange.at(ncp);
      auto rMin = rR.lower.value_or(0.);
      auto rMax = rR.upper.value_or(std::numeric_limits<double>::max());
      auto zR = volume.zRange.at(ncp);
      auto zMin = zR.lower.value_or(-std::numeric_limits<double>::max());
      auto zMax = zR.upper.value_or(std::numeric_limits<double>::max());
      lConfig.parseRanges = {
          {AxisDirection::AxisR, {rMin, rMax}},
          {AxisDirection::AxisZ, {zMin, zMax}},
      };

      // Fill the layer splitting parameters in r/z
      auto str = volume.splitTolR.at(ncp);
      auto stz = volume.splitTolZ.at(ncp);
      if (0 < str) {
        lConfig.splitConfigs.emplace_back(AxisDirection::AxisR, str);
      }
      if (0 < stz) {
        lConfig.splitConfigs.emplace_back(AxisDirection::AxisZ, stz);
      }
      lConfig.binning0 = volume.binning0.at(ncp);
      lConfig.binning1 = volume.binning1.at(ncp);

      layerBuilderConfig.layerConfigurations[ncp].push_back(lConfig);
    }

    // Perform splitting of cylinders and discs
    if (volume.cylinderDiscSplit) {
      TGeoCylinderDiscSplitter::Config cdsConfig;
      cdsConfig.cylinderPhiSegments = volume.cylinderNPhiSegments;
      cdsConfig.cylinderLongitudinalSegments = volume.cylinderNZSegments;
      cdsConfig.discPhiSegments = volume.discNPhiSegments;
      cdsConfig.discRadialSegments = volume.discNRSegments;
      layerBuilderConfig.detectorElementSplitter =
          std::make_shared<const TGeoCylinderDiscSplitter>(
              cdsConfig,
              logger.clone("TGeoCylinderDiscSplitter", config.layerLogLevel));
    } else if (volume.itkModuleSplit) {
      TGeoITkModuleSplitter::Config itkConfig;
      itkConfig.barrelMap = volume.barrelMap;
      itkConfig.discMap = volume.discMap;
      itkConfig.splitPatterns = volume.splitPatterns;
      layerBuilderConfig.detectorElementSplitter =
          std::make_shared<TGeoITkModuleSplitter>(
              itkConfig,
              logger.clone("TGeoITkModuleSplitter", config.layerLogLevel));
    }

    detLayerConfigs.push_back(layerBuilderConfig);
  }

  return detLayerConfigs;
}

/// @brief Function to build the generic tracking geometry
// from a TGeo object.
///
/// It does *currently* not translate the material, this has
/// to be done with a material mapping stage
///
/// @tparam variable_map_t is the variable map
///
/// @param vm is the variable map from the options
std::shared_ptr<const TrackingGeometry> buildTGeoDetector(
    const TGeoDetector::Config& config, const GeometryContext& context,
    std::vector<std::shared_ptr<const DetectorElementBase>>& detElementStore,
    std::shared_ptr<const IMaterialDecorator> materialDecorator,
    const Logger& logger) {
  // configure surface array creator
  SurfaceArrayCreator::Config sacConfig;
  auto surfaceArrayCreator = std::make_shared<const SurfaceArrayCreator>(
      sacConfig, logger.clone("SurfaceArrayCreator", config.surfaceLogLevel));
  // configure the proto layer helper
  ProtoLayerHelper::Config plhConfig;
  auto protoLayerHelper = std::make_shared<const ProtoLayerHelper>(
      plhConfig, logger.clone("ProtoLayerHelper", config.layerLogLevel));
  // configure the layer creator that uses the surface array creator
  LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<const LayerCreator>(
      lcConfig, logger.clone("LayerCreator", config.layerLogLevel));
  // configure the layer array creator
  LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const LayerArrayCreator>(
      lacConfig, logger.clone("LayerArrayCreator", config.layerLogLevel));
  // tracking volume array creator
  TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator = std::make_shared<const TrackingVolumeArrayCreator>(
      tvacConfig,
      logger.clone("TrackingVolumeArrayCreator", config.volumeLogLevel));
  // configure the cylinder volume helper
  CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper = std::make_shared<const CylinderVolumeHelper>(
      cvhConfig, logger.clone("CylinderVolumeHelper", config.volumeLogLevel));

  //-------------------------------------------------------------------------------------
  // list the volume builders
  std::list<std::shared_ptr<const ITrackingVolumeBuilder>> volumeBuilders;

  // Create a beam pipe if configured to do so
  if (config.buildBeamPipe) {
    /// configure the beam pipe layer builder
    PassiveLayerBuilder::Config bplConfig;
    bplConfig.layerIdentification = "BeamPipe";
    bplConfig.centralLayerRadii = {config.beamPipeRadius};
    bplConfig.centralLayerHalflengthZ = {config.beamPipeHalflengthZ};
    bplConfig.centralLayerThickness = {config.beamPipeLayerThickness};
    auto beamPipeBuilder = std::make_shared<const PassiveLayerBuilder>(
        bplConfig, logger.clone("BeamPipeLayerBuilder", config.layerLogLevel));
    // create the volume for the beam pipe
    CylinderVolumeBuilder::Config bpvConfig;
    bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
    bpvConfig.volumeName = "BeamPipe";
    bpvConfig.layerBuilder = beamPipeBuilder;
    bpvConfig.layerEnvelopeR = {config.beamPipeEnvelopeR,
                                config.beamPipeEnvelopeR};
    bpvConfig.buildToRadiusZero = true;
    auto beamPipeVolumeBuilder = std::make_shared<const CylinderVolumeBuilder>(
        bpvConfig,
        logger.clone("BeamPipeVolumeBuilder", config.volumeLogLevel));
    // add to the list of builders
    volumeBuilders.push_back(beamPipeVolumeBuilder);
  }

  // Import the file from
  TGeoManager::Import(config.fileName.c_str());

  auto layerBuilderConfigs = makeLayerBuilderConfigs(config, logger);

  // Remember the layer builders to collect the detector elements
  std::vector<std::shared_ptr<const ActsPlugins::TGeoLayerBuilder>>
      tgLayerBuilders;

  for (auto& lbc : layerBuilderConfigs) {
    std::shared_ptr<const LayerCreator> layerCreatorLB = nullptr;

    if (lbc.autoSurfaceBinning) {
      // Configure surface array creator (optionally) per layer builder
      // (in order to configure them to work appropriately)
      SurfaceArrayCreator::Config sacConfigLB;
      sacConfigLB.surfaceMatcher = lbc.surfaceBinMatcher;
      auto surfaceArrayCreatorLB = std::make_shared<const SurfaceArrayCreator>(
          sacConfigLB,
          logger.clone(lbc.configurationName + "SurfaceArrayCreator",
                       config.surfaceLogLevel));
      // configure the layer creator that uses the surface array creator
      LayerCreator::Config lcConfigLB;
      lcConfigLB.surfaceArrayCreator = surfaceArrayCreatorLB;
      layerCreatorLB = std::make_shared<const LayerCreator>(
          lcConfigLB, logger.clone(lbc.configurationName + "LayerCreator",
                                   config.layerLogLevel));
    }

    // Configure the proto layer helper
    ProtoLayerHelper::Config plhConfigLB;
    auto protoLayerHelperLB = std::make_shared<const ProtoLayerHelper>(
        plhConfigLB, logger.clone(lbc.configurationName + "ProtoLayerHelper",
                                  config.layerLogLevel));

    //-------------------------------------------------------------------------------------
    lbc.layerCreator =
        (layerCreatorLB != nullptr) ? layerCreatorLB : layerCreator;
    lbc.protoLayerHelper =
        (protoLayerHelperLB != nullptr) ? protoLayerHelperLB : protoLayerHelper;

    auto layerBuilder = std::make_shared<const ActsPlugins::TGeoLayerBuilder>(
        lbc, logger.clone(lbc.configurationName + "LayerBuilder",
                          config.layerLogLevel));
    // remember the layer builder
    tgLayerBuilders.push_back(layerBuilder);

    // build the pixel volume
    CylinderVolumeBuilder::Config volumeConfig;
    volumeConfig.trackingVolumeHelper = cylinderVolumeHelper;
    volumeConfig.volumeName = lbc.configurationName;
    volumeConfig.buildToRadiusZero = volumeBuilders.empty();
    volumeConfig.layerEnvelopeR = {config.layerEnvelopeR,
                                   config.layerEnvelopeR};
    auto ringLayoutConfiguration =
        [&](const std::vector<ActsPlugins::TGeoLayerBuilder::LayerConfig>&
                lConfigs) {
          for (const auto& lcfg : lConfigs) {
            for (const auto& scfg : lcfg.splitConfigs) {
              if (scfg.first == AxisDirection::AxisR && scfg.second > 0.) {
                volumeConfig.ringTolerance =
                    std::max(volumeConfig.ringTolerance, scfg.second);
                volumeConfig.checkRingLayout = true;
              }
            }
          }
        };
    ringLayoutConfiguration(lbc.layerConfigurations[0]);
    ringLayoutConfiguration(lbc.layerConfigurations[2]);
    volumeConfig.layerBuilder = layerBuilder;
    auto volumeBuilder = std::make_shared<const CylinderVolumeBuilder>(
        volumeConfig, logger.clone(lbc.configurationName + "VolumeBuilder",
                                   config.volumeLogLevel));
    // add to the list of builders
    volumeBuilders.push_back(volumeBuilder);
  }

  //-------------------------------------------------------------------------------------
  // create the tracking geometry
  TrackingGeometryBuilder::Config tgConfig;
  // Add the builders
  tgConfig.materialDecorator = std::move(materialDecorator);
  tgConfig.geometryIdentifierHook = config.geometryIdentifierHook;

  for (auto& vb : volumeBuilders) {
    tgConfig.trackingVolumeBuilders.push_back(
        [=](const auto& gcontext, const auto& inner, const auto&) {
          return vb->trackingVolume(gcontext, inner);
        });
  }
  // Add the helper
  tgConfig.trackingVolumeHelper = cylinderVolumeHelper;
  auto cylinderGeometryBuilder =
      std::make_shared<const TrackingGeometryBuilder>(
          tgConfig,
          logger.clone("TrackerGeometryBuilder", config.volumeLogLevel));
  // get the geometry
  auto trackingGeometry = cylinderGeometryBuilder->trackingGeometry(context);
  // collect the detector element store
  for (auto& lBuilder : tgLayerBuilders) {
    auto detElements = lBuilder->detectorElements();
    detElementStore.insert(detElementStore.begin(), detElements.begin(),
                           detElements.end());
  }

  /// return the tracking geometry
  return trackingGeometry;
}

}  // namespace

/// Read the TGeo layer builder configurations from the user configuration.
void TGeoDetector::readTGeoLayerBuilderConfigsFile(const std::string& path,
                                                   Config& config) {
  if (path.empty()) {
    return;
  }
  nlohmann::json djson;
  std::ifstream infile(path, std::ifstream::in | std::ifstream::binary);
  infile >> djson;

  config.unitScalor = djson["geo-tgeo-unit-scalor"];

  config.buildBeamPipe = djson["geo-tgeo-build-beampipe"];
  if (config.buildBeamPipe) {
    const auto beamPipeParameters =
        djson["geo-tgeo-beampipe-parameters"].get<std::array<double, 3>>();
    config.beamPipeRadius = beamPipeParameters[0];
    config.beamPipeHalflengthZ = beamPipeParameters[1];
    config.beamPipeLayerThickness = beamPipeParameters[2];
  }

  // Fill nested volume configs
  for (const auto& volume : djson["Volumes"]) {
    auto& vol = config.volumes.emplace_back();
    vol = volume;
  }
}

TGeoDetector::TGeoDetector(const Config& cfg)
    : Detector(getDefaultLogger("TGeoDetector", cfg.logLevel)), m_cfg(cfg) {
  m_nominalGeometryContext = GeometryContext::dangerouslyDefaultConstruct();

  m_trackingGeometry =
      buildTGeoDetector(m_cfg, m_nominalGeometryContext, m_detectorStore,
                        m_cfg.materialDecorator, logger());
}

void TGeoDetector::Config::readJson(const std::string& jsonFile) {
  readTGeoLayerBuilderConfigsFile(jsonFile, *this);
}
std::shared_ptr<const TrackingGeometry> buildTGeoDetectorWrapper(
    const TGeoDetector::Config& config, const GeometryContext& context,
    std::vector<std::shared_ptr<const DetectorElementBase>>& detElementStore,
    std::shared_ptr<const IMaterialDecorator> materialDecorator,
    const Logger& logger) {
  return buildTGeoDetector(config, context, detElementStore,
                           std::move(materialDecorator), logger);
}
}  // namespace ActsExamples
