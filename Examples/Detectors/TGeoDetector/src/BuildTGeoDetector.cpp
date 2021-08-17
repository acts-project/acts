// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TGeoDetector/BuildTGeoDetector.hpp"

#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/SurfaceBinningMatcher.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Plugins/TGeo/TGeoCylinderDiscSplitter.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsExamples/TGeoDetector/BuildTGeoDetector.hpp"
#include "ActsExamples/TGeoDetector/TGeoDetectorOptions.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <list>
#include <vector>

#include <TGeoManager.h>

/// @brief global method to build the generic tracking geometry
// from a TGeo object.
///
/// It does *currently* not translate the material, this has
/// to be done with a material mapping stage
///
/// @tparam variable_map_t is the variable map
///
/// @param vm is the variable map from the options
std::shared_ptr<const Acts::TrackingGeometry>
ActsExamples::TGeo::buildTGeoDetector(
    const TGeoDetector::Config& config, const Acts::GeometryContext& context,
    std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>&
        detElementStore,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) {
  // configure surface array creator
  Acts::SurfaceArrayCreator::Config sacConfig;
  auto surfaceArrayCreator = std::make_shared<const Acts::SurfaceArrayCreator>(
      sacConfig,
      Acts::getDefaultLogger("SurfaceArrayCreator", config.surfaceLogLevel));
  // configure the proto layer helper
  Acts::ProtoLayerHelper::Config plhConfig;
  auto protoLayerHelper = std::make_shared<const Acts::ProtoLayerHelper>(
      plhConfig,
      Acts::getDefaultLogger("ProtoLayerHelper", config.layerLogLevel));
  // configure the layer creator that uses the surface array creator
  Acts::LayerCreator::Config lcConfig;
  lcConfig.surfaceArrayCreator = surfaceArrayCreator;
  auto layerCreator = std::make_shared<const Acts::LayerCreator>(
      lcConfig, Acts::getDefaultLogger("LayerCreator", config.layerLogLevel));
  // configure the layer array creator
  Acts::LayerArrayCreator::Config lacConfig;
  auto layerArrayCreator = std::make_shared<const Acts::LayerArrayCreator>(
      lacConfig,
      Acts::getDefaultLogger("LayerArrayCreator", config.layerLogLevel));
  // tracking volume array creator
  Acts::TrackingVolumeArrayCreator::Config tvacConfig;
  auto tVolumeArrayCreator =
      std::make_shared<const Acts::TrackingVolumeArrayCreator>(
          tvacConfig, Acts::getDefaultLogger("TrackingVolumeArrayCreator",
                                             config.volumeLogLevel));
  // configure the cylinder volume helper
  Acts::CylinderVolumeHelper::Config cvhConfig;
  cvhConfig.layerArrayCreator = layerArrayCreator;
  cvhConfig.trackingVolumeArrayCreator = tVolumeArrayCreator;
  auto cylinderVolumeHelper =
      std::make_shared<const Acts::CylinderVolumeHelper>(
          cvhConfig, Acts::getDefaultLogger("CylinderVolumeHelper",
                                            config.volumeLogLevel));

  //-------------------------------------------------------------------------------------
  // list the volume builders
  std::list<std::shared_ptr<const Acts::ITrackingVolumeBuilder>> volumeBuilders;

  // Create a beam pipe if configured to do so
  if (config.buildBeamPipe) {
    /// configure the beam pipe layer builder
    Acts::PassiveLayerBuilder::Config bplConfig;
    bplConfig.layerIdentification = "BeamPipe";
    bplConfig.centralLayerRadii = std::vector<double>(1, config.beamPipeRadius);
    bplConfig.centralLayerHalflengthZ =
        std::vector<double>(1, config.beamPipeHalflengthZ);
    bplConfig.centralLayerThickness =
        std::vector<double>(1, config.beamPipeLayerThickness);
    auto beamPipeBuilder = std::make_shared<const Acts::PassiveLayerBuilder>(
        bplConfig,
        Acts::getDefaultLogger("BeamPipeLayerBuilder", config.layerLogLevel));
    // create the volume for the beam pipe
    Acts::CylinderVolumeBuilder::Config bpvConfig;
    bpvConfig.trackingVolumeHelper = cylinderVolumeHelper;
    bpvConfig.volumeName = "BeamPipe";
    bpvConfig.layerBuilder = beamPipeBuilder;
    bpvConfig.layerEnvelopeR = {1. * Acts::UnitConstants::mm,
                                1. * Acts::UnitConstants::mm};
    bpvConfig.buildToRadiusZero = true;
    auto beamPipeVolumeBuilder =
        std::make_shared<const Acts::CylinderVolumeBuilder>(
            bpvConfig, Acts::getDefaultLogger("BeamPipeVolumeBuilder",
                                              config.volumeLogLevel));
    // add to the list of builders
    volumeBuilders.push_back(beamPipeVolumeBuilder);
  }

  // Import the file from
  TGeoManager::Import(config.fileName.c_str());

  auto layerBuilderConfigs =
      ActsExamples::TGeo::makeLayerBuilderConfigs(config);

  // Remember the layer builders to collect the detector elements
  std::vector<std::shared_ptr<const Acts::TGeoLayerBuilder>> tgLayerBuilders;

  for (auto& lbc : layerBuilderConfigs) {
    std::shared_ptr<const Acts::LayerCreator> layerCreatorLB = nullptr;

    if (lbc.autoSurfaceBinning) {
      // Configure surface array creator (optionally) per layer builder
      // (in order to configure them to work appropriately)
      Acts::SurfaceArrayCreator::Config sacConfigLB;
      sacConfigLB.surfaceMatcher = lbc.surfaceBinMatcher;
      auto surfaceArrayCreatorLB =
          std::make_shared<const Acts::SurfaceArrayCreator>(
              sacConfigLB, Acts::getDefaultLogger(
                               lbc.configurationName + "SurfaceArrayCreator",
                               config.surfaceLogLevel));
      // configure the layer creator that uses the surface array creator
      Acts::LayerCreator::Config lcConfigLB;
      lcConfigLB.surfaceArrayCreator = surfaceArrayCreatorLB;
      layerCreatorLB = std::make_shared<const Acts::LayerCreator>(
          lcConfigLB,
          Acts::getDefaultLogger(lbc.configurationName + "LayerCreator",
                                 config.layerLogLevel));
    }

    // Configure the proto layer helper
    Acts::ProtoLayerHelper::Config plhConfigLB;
    auto protoLayerHelperLB = std::make_shared<const Acts::ProtoLayerHelper>(
        plhConfigLB,
        Acts::getDefaultLogger(lbc.configurationName + "ProtoLayerHelper",
                               config.layerLogLevel));

    //-------------------------------------------------------------------------------------
    lbc.layerCreator =
        (layerCreatorLB != nullptr) ? layerCreatorLB : layerCreator;
    lbc.protoLayerHelper =
        (protoLayerHelperLB != nullptr) ? protoLayerHelperLB : protoLayerHelper;

    auto layerBuilder = std::make_shared<const Acts::TGeoLayerBuilder>(
        lbc, Acts::getDefaultLogger(lbc.configurationName + "LayerBuilder",
                                    config.layerLogLevel));
    // remember the layer builder
    tgLayerBuilders.push_back(layerBuilder);

    // build the pixel volume
    Acts::CylinderVolumeBuilder::Config volumeConfig;
    volumeConfig.trackingVolumeHelper = cylinderVolumeHelper;
    volumeConfig.volumeName = lbc.configurationName;
    volumeConfig.buildToRadiusZero = (volumeBuilders.size() == 0);
    volumeConfig.layerEnvelopeR = {1. * Acts::UnitConstants::mm,
                                   5. * Acts::UnitConstants::mm};
    auto ringLayoutConfiguration =
        [&](const std::vector<Acts::TGeoLayerBuilder::LayerConfig>& lConfigs)
        -> void {
      for (const auto& lcfg : lConfigs) {
        for (const auto& scfg : lcfg.splitConfigs) {
          if (scfg.first == Acts::binR and scfg.second > 0.) {
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
    volumeConfig.volumeSignature = 0;
    auto volumeBuilder = std::make_shared<const Acts::CylinderVolumeBuilder>(
        volumeConfig,
        Acts::getDefaultLogger(lbc.configurationName + "VolumeBuilder",
                               config.volumeLogLevel));
    // add to the list of builders
    volumeBuilders.push_back(volumeBuilder);
  }

  //-------------------------------------------------------------------------------------
  // create the tracking geometry
  Acts::TrackingGeometryBuilder::Config tgConfig;
  // Add the builders
  tgConfig.materialDecorator = mdecorator;

  for (auto& vb : volumeBuilders) {
    tgConfig.trackingVolumeBuilders.push_back(
        [=](const auto& gcontext, const auto& inner, const auto&) {
          return vb->trackingVolume(gcontext, inner);
        });
  }
  // Add the helper
  tgConfig.trackingVolumeHelper = cylinderVolumeHelper;
  auto cylinderGeometryBuilder =
      std::make_shared<const Acts::TrackingGeometryBuilder>(
          tgConfig, Acts::getDefaultLogger("TrackerGeometryBuilder",
                                           config.volumeLogLevel));
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

std::vector<Acts::TGeoLayerBuilder::Config>
ActsExamples::TGeo::makeLayerBuilderConfigs(
    const TGeoDetector::Config& config) {
  std::vector<Acts::TGeoLayerBuilder::Config> detLayerConfigs;

  // iterate over all configured detector volumes
  for (const auto& volume : config.volumes) {
    Acts::TGeoLayerBuilder::Config layerBuilderConfig;
    layerBuilderConfig.configurationName = volume.name;
    layerBuilderConfig.unit = config.unitScalor;

    // configure surface autobinning
    std::vector<std::pair<double, double>> binTolerances(
        static_cast<size_t>(Acts::binValues), {0., 0.});
    binTolerances[Acts::binR] = {volume.binToleranceR.lower.value_or(0.),
                                 volume.binToleranceR.upper.value_or(0.)};
    binTolerances[Acts::binZ] = {volume.binToleranceZ.lower.value_or(0.),
                                 volume.binToleranceZ.upper.value_or(0.)};
    binTolerances[Acts::binPhi] = {volume.binTolerancePhi.lower.value_or(0.),
                                   volume.binTolerancePhi.upper.value_or(0.)};
    layerBuilderConfig.autoSurfaceBinning = true;
    layerBuilderConfig.surfaceBinMatcher =
        Acts::SurfaceBinningMatcher(binTolerances);

    // loop over the negative/central/positive layer configurations
    for (auto ncp : {
             TGeoDetector::Config::Negative,
             TGeoDetector::Config::Central,
             TGeoDetector::Config::Positive,
         }) {
      if (!volume.layers.at(ncp)) {
        continue;
      }

      Acts::TGeoLayerBuilder::LayerConfig lConfig;
      lConfig.volumeName = volume.subVolumeName.at(ncp);
      lConfig.sensorNames = volume.sensitiveNames.at(ncp);
      lConfig.localAxes = volume.sensitiveAxes.at(ncp);

      // Fill the parsing restrictions in r/z
      auto rR = volume.rRange.at(ncp);
      auto rMin = rR.lower.value_or(0.);
      auto rMax = rR.upper.value_or(std::numeric_limits<double>::max());
      auto zR = volume.zRange.at(ncp);
      auto zMin = zR.lower.value_or(-std::numeric_limits<double>::max());
      auto zMax = zR.upper.value_or(std::numeric_limits<double>::max());
      lConfig.parseRanges = {
          {Acts::binR, {rMin, rMax}},
          {Acts::binZ, {zMin, zMax}},
      };

      // Fill the layer splitting parameters in r/z
      auto str = volume.splitTolR.at(ncp);
      auto stz = volume.splitTolZ.at(ncp);
      if (0 < str) {
        lConfig.splitConfigs.emplace_back(Acts::binR, str);
      }
      if (0 < stz) {
        lConfig.splitConfigs.emplace_back(Acts::binZ, stz);
      }

      layerBuilderConfig.layerConfigurations[ncp].push_back(lConfig);
    }

    // Perform splitting of cylinders and discs
    if (volume.cylinderDiscSplit) {
      Acts::TGeoCylinderDiscSplitter::Config cdsConfig;
      cdsConfig.cylinderPhiSegments = volume.cylinderNPhiSegments;
      cdsConfig.cylinderLongitudinalSegments = volume.cylinderNZSegments;
      cdsConfig.discPhiSegments = volume.discNPhiSegments;
      cdsConfig.discRadialSegments = volume.discNRSegments;
      layerBuilderConfig.detectorElementSplitter =
          std::make_shared<const Acts::TGeoCylinderDiscSplitter>(cdsConfig);
    }

    detLayerConfigs.push_back(layerBuilderConfig);
  }

  return detLayerConfigs;
}
