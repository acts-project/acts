// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdio.h>

#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include "Acts/Plugins/TGeo/TGeoParser.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelper.hpp"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

Acts::TGeoLayerBuilder::TGeoLayerBuilder(
    const Acts::TGeoLayerBuilder::Config& config,
    std::unique_ptr<const Logger> logger)
    : m_cfg(), m_logger(std::move(logger)) {
  setConfiguration(config);
}

Acts::TGeoLayerBuilder::~TGeoLayerBuilder() = default;

void Acts::TGeoLayerBuilder::setConfiguration(
    const Acts::TGeoLayerBuilder::Config& config) {
  m_cfg = config;
}

void Acts::TGeoLayerBuilder::setLogger(
    std::unique_ptr<const Logger> newLogger) {
  m_logger = std::move(newLogger);
}

const Acts::LayerVector Acts::TGeoLayerBuilder::negativeLayers(
    const GeometryContext& gctx) const {
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector nVector;
  mutableThis->buildLayers(gctx, nVector, -1);
  return nVector;
}

const Acts::LayerVector Acts::TGeoLayerBuilder::centralLayers(
    const GeometryContext& gctx) const {
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector cVector;
  mutableThis->buildLayers(gctx, cVector, 0);
  return cVector;
}

const Acts::LayerVector Acts::TGeoLayerBuilder::positiveLayers(
    const GeometryContext& gctx) const {
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector pVector;
  mutableThis->buildLayers(gctx, pVector, 1);
  return pVector;
}

void Acts::TGeoLayerBuilder::buildLayers(const GeometryContext& gctx,
                                         LayerVector& layers, int type) {
  // Bail out if you have no gGeoManager
  if (gGeoManager == nullptr) {
    ACTS_WARNING("No gGeoManager found - bailing out.");
    return;
  }

  using LayerSurfaceVector = std::vector<std::shared_ptr<const Surface>>;
  LayerSurfaceVector layerSurfaces;

  std::vector<LayerConfig> layerConfigs = m_cfg.layerConfigurations[type + 1];
  std::string layerType = m_layerTypes[type + 1];

  // Appropriate screen output
  std::string addonOutput = m_cfg.layerSplitToleranceR[type + 1] > 0.
                                ? std::string(", splitting in r")
                                : std::string("");
  addonOutput += m_cfg.layerSplitToleranceZ[type + 1] > 0.
                     ? std::string(", splitting in z")
                     : std::string("");
  addonOutput += std::string(".");

  // Screen output of the configuration
  ACTS_DEBUG(layerType << " layers : found " << layerConfigs.size()
                       << " configuration(s)" + addonOutput);

  // Helper function to fill the layer
  auto fillLayer = [&](const LayerSurfaceVector lSurfaces,
                       const LayerConfig& lCfg) -> void {
    // Create the layer  - either way as cylinder or disk
    if (type == 0) {
      ACTS_DEBUG("- creating CylinderLayer with " << lSurfaces.size()
                                                  << " surfaces.");
      ProtoLayer pl(gctx, lSurfaces);
      pl.envelope[Acts::binR] = {lCfg.envelope.first, lCfg.envelope.second};
      pl.envelope[Acts::binZ] = {lCfg.envelope.second, lCfg.envelope.second};
      layers.push_back(m_cfg.layerCreator->cylinderLayer(
          gctx, lSurfaces, lCfg.binsLoc0, lCfg.binsLoc1, pl));
    } else {
      ACTS_DEBUG("- creating DiscLayer with " << lSurfaces.size()
                                              << " surfaces.");
      ProtoLayer pl(gctx, lSurfaces);
      pl.envelope[Acts::binR] = {lCfg.envelope.first, lCfg.envelope.second};
      pl.envelope[Acts::binZ] = {lCfg.envelope.second, lCfg.envelope.second};
      layers.push_back(m_cfg.layerCreator->discLayer(
          gctx, lSurfaces, lCfg.binsLoc0, lCfg.binsLoc1, pl));
    }
  };

  for (auto layerCfg : layerConfigs) {
    ACTS_DEBUG("- layer configuration found for layer " << layerCfg.layerName
                                                        << " with sensors ");
    for (auto& sensor : layerCfg.sensorNames) {
      ACTS_DEBUG("  - sensor: " << sensor);
    }
    if (not layerCfg.parseRanges.empty()) {
      for (const auto& pRange : layerCfg.parseRanges) {
        ACTS_DEBUG("- layer parsing restricted in "
                   << binningValueNames[pRange.first] << " to ["
                   << pRange.second.first << "/" << pRange.second.second
                   << "].");
      }
    }
    if (not layerCfg.splitConfigs.empty()) {
      for (const auto& sConfig : layerCfg.splitConfigs) {
        ACTS_DEBUG("- layer splitting attempt in "
                   << binningValueNames[sConfig.first] << " with tolerance "
                   << sConfig.second << ".");
      }
    }

    // Step down from the top volume each time to collect the logical tree
    TGeoVolume* tvolume = gGeoManager->GetTopVolume();
    if (tvolume != nullptr) {
      TGeoParser::Options tgpOptions;
      tgpOptions.volumeNames = {layerCfg.layerName};
      tgpOptions.targetNames = layerCfg.sensorNames;
      tgpOptions.parseRanges = layerCfg.parseRanges;
      tgpOptions.unit = m_cfg.unit;

      TGeoParser::State tgpState;
      tgpState.volume = tvolume;

      TGeoParser::select(tgpState, tgpOptions);

      ACTS_DEBUG("- number of selsected nodes found : "
                 << tgpState.selectedNodes.size());

      for (auto& snode : tgpState.selectedNodes) {
        auto identifier =
            m_cfg.identifierProvider != nullptr
                ? m_cfg.identifierProvider->identify(gctx, *snode.node)
                : Identifier();

        auto tgElement = std::make_shared<const Acts::TGeoDetectorElement>(
            identifier, *snode.node, *snode.transform, layerCfg.localAxes,
            m_cfg.unit);
        m_elementStore.push_back(tgElement);
        layerSurfaces.push_back(tgElement->surface().getSharedPtr());
      }

      ACTS_DEBUG("- created TGeoDetectorElements : " << layerSurfaces.size());

      if (m_cfg.protoLayerHelper != nullptr and
          not layerCfg.splitConfigs.empty()) {
        auto protoLayers = m_cfg.protoLayerHelper->protoLayers(
            gctx, unpack_shared_vector(layerSurfaces), layerCfg.splitConfigs);
        ACTS_DEBUG("- splitting into " << protoLayers.size() << " layers.");
        for (auto& pLayer : protoLayers) {
          layerSurfaces.clear();
          for (const auto& lsurface : pLayer.surfaces()) {
            layerSurfaces.push_back(lsurface->getSharedPtr());
          }
          fillLayer(layerSurfaces, layerCfg);
        }
      } else {
        fillLayer(layerSurfaces, layerCfg);
      }
    }
  }
  return;
}
