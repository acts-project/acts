// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include <stdio.h>
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
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
  mutableThis->buildLayers(gctx, pVector, -1);
  return pVector;
}

void Acts::TGeoLayerBuilder::buildLayers(const GeometryContext& gctx,
                                         LayerVector& layers, int type) {
  // bail out if you have no gGeoManager
  if (gGeoManager == nullptr) {
    return;
  }

  // Prepare which ones to build
  std::vector<LayerConfig> layerConfigs;
  std::string layerType = "No";
  switch (type) {
    case -1: {
      layerConfigs = m_cfg.negativeLayerConfigs;
      layerType = "Negative";
    } break;
    case 0: {
      layerConfigs = m_cfg.centralLayerConfigs;
      layerType = "Central";
    } break;
    case 1: {
      layerConfigs = m_cfg.positiveLayerConfigs;
      layerType = "Positive";
    } break;
  }
  // screen output
  ACTS_DEBUG(layerType << " Layers : found " << layerConfigs.size()
                       << " configurations.");
  for (auto layerCfg : layerConfigs) {
    // prepare the layer surfaces
    using LayerSurfaceVector = std::vector<std::shared_ptr<const Surface>>;
    LayerSurfaceVector layerSurfaces;

    ACTS_DEBUG("- layer configuration found for layer "
               << layerCfg.layerName << " with sensor " << layerCfg.sensorName);
    // we have to step down from the top volume each time to collect the logical
    // tree
    TGeoVolume* tvolume = gGeoManager->GetTopVolume();
    if (tvolume != nullptr) {
      // recursively step down
      resolveSensitive(gctx, layerSurfaces, tvolume, nullptr, TGeoIdentity(),
                       layerCfg, type);
      // screen output
      ACTS_DEBUG(
          "- number of senstive sensors found : " << layerSurfaces.size());

      // Helper function to fill the layer
      auto fillLayer = [&](const LayerSurfaceVector lSurfaces,
                           const LayerConfig& lCfg) -> void {
        // create the layer  - either way
        if (type == 0) {
          ProtoLayer pl(gctx, lSurfaces);
          pl.envR = {lCfg.envelope.first, lCfg.envelope.second};
          pl.envZ = {lCfg.envelope.second, lCfg.envelope.second};
          layers.push_back(m_cfg.layerCreator->cylinderLayer(
              gctx, lSurfaces, lCfg.binsLoc0, lCfg.binsLoc1, pl));
        } else {
          ProtoLayer pl(gctx, lSurfaces);
          pl.envR = {lCfg.envelope.first, lCfg.envelope.second};
          pl.envZ = {lCfg.envelope.second, lCfg.envelope.second};
          layers.push_back(m_cfg.layerCreator->discLayer(
              gctx, lSurfaces, lCfg.binsLoc0, lCfg.binsLoc1, pl));
        }
      };

      // Check if a radial split is requested
      if (layerCfg.splitRadii.size() > 1) {
        ACTS_DEBUG("- radially split layers seperated by more than "
                   << m_cfg.centralLayerSplit);
        ACTS_DEBUG("- surface center r min/max = " << layerCfg.rminmax.first
                                                   << ", "
                                                   << layerCfg.rminmax.second);
        ACTS_DEBUG("- number of proposed split radii is "
                   << layerCfg.splitRadii.size());

        // Prepare the vector of split surfaces
        std::vector<LayerSurfaceVector> splitSurfaces{
            layerCfg.splitRadii.size(), LayerSurfaceVector{}};
        for (const auto& surface : layerSurfaces) {
          double surfaceR = surface->binningPositionValue(gctx, binR);
          unsigned ir = 0;
          for (const auto& sugr : layerCfg.splitRadii) {
            if (std::abs(sugr - surfaceR) < m_cfg.centralLayerSplit) {
              splitSurfaces[ir].push_back(surface);
            }
            ++ir;
          }
        }

        ACTS_DEBUG("Result of the split analysis:");
        unsigned il = 0;
        for (const auto& lSurfaces : splitSurfaces) {
          ACTS_DEBUG("- layer  " << il << " has " << lSurfaces.size()
                                 << " surfaces.");
          fillLayer(lSurfaces, layerCfg);
        }
        return;
      }
      // No splitting done
      fillLayer(layerSurfaces, layerCfg);
    }
  }
}

void Acts::TGeoLayerBuilder::resolveSensitive(
    const GeometryContext& gctx,
    std::vector<std::shared_ptr<const Acts::Surface>>& layerSurfaces,
    TGeoVolume* tgVolume, TGeoNode* tgNode, const TGeoMatrix& tgTransform,
    LayerConfig& layerConfig, int type, bool correctBranch,
    const std::string& offset) {
  /// some screen output for disk debugging
  if (type != 0) {
    const Double_t* ctranslation = tgTransform.GetTranslation();
    ACTS_VERBOSE(offset << "current z translation is : " << ctranslation[2]);
  }

  if (tgVolume != nullptr) {
    std::string volumeName = tgVolume->GetName();
    /// some screen output indicating that the volume was found
    ACTS_VERBOSE(offset << "[o] Volume : " << volumeName
                        << " - checking for volume name "
                        << layerConfig.layerName);

    // once in the current branch, always in the current branch
    bool correctVolume = correctBranch;
    if (!correctVolume &&
        (volumeName.find(layerConfig.layerName) != std::string::npos ||
         match(layerConfig.layerName.c_str(), volumeName.c_str()))) {
      correctVolume = true;
      ACTS_VERBOSE(offset << "    triggered current branch!");
    }
    // loop over the daughters and collect them
    auto daugthers = tgVolume->GetNodes();
    // screen output
    ACTS_VERBOSE(offset << "has " << tgVolume->GetNdaughters()
                        << " daughters.");
    // a daughter iterator
    TIter iObj(daugthers);
    // while loop over the objects
    while (TObject* obj = iObj()) {
      // dynamic_cast to a node
      TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
      if (node != nullptr) {
        resolveSensitive(gctx, layerSurfaces, nullptr, node, tgTransform,
                         layerConfig, type, correctVolume, offset + "  ");
      }
    }
  } else {
    ACTS_VERBOSE("No volume present.");
  }

  /// if you have a node, get the volume and step down further
  if (tgNode != nullptr) {
    // get the matrix of the current
    const TGeoMatrix* tgMatrix = tgNode->GetMatrix();
    /// get the translation of the parent
    const Double_t* translation = tgTransform.GetTranslation();
    // get the z value
    double z = translation[2];
    // get the name of the node
    std::string tNodeName = tgNode->GetName();
    ACTS_VERBOSE(offset << "[>] Node : " << tNodeName
                        << " - checking for sensor name "
                        << layerConfig.sensorName);
    // find out the branch hit - single layer depth is supported by
    // sensor==layer
    bool branchHit =
        correctBranch || (layerConfig.sensorName == layerConfig.layerName);
    if (branchHit &&
        (tNodeName.find(layerConfig.sensorName) != std::string::npos ||
         match(layerConfig.sensorName.c_str(), tNodeName.c_str()))) {
      ACTS_VERBOSE(offset << "Sensor name found in correct branch.");

      // set the visibility to kTrue
      if (m_cfg.setVisibility) {
        tgNode->SetVisibility(kTRUE);
      }
      // create the detector element - check on the type for the size
      if ((type == 0) || type * z > 0.) {
        //  senstive volume found, collect it
        ACTS_VERBOSE(offset << "[>>] accepted !");
        // create the element
        auto tgElement = std::make_shared<const Acts::TGeoDetectorElement>(
            Identifier(), tgNode, &tgTransform, layerConfig.localAxes,
            m_cfg.unit);
        // record the element @todo solve with provided cache
        m_elementStore.push_back(tgElement);
        // record the surface
        // element owns the surface, we give shared ownership to
        // layer surfaces -> i.e. the Layer to be built
        layerSurfaces.push_back(tgElement->surface().getSharedPtr());
        // Record rmin/rmax in the layerConfig for eventual splitting
        double surfaceR = tgElement->surface().binningPositionValue(gctx, binR);
        layerConfig.rminmax.first =
            std::min(layerConfig.rminmax.first, surfaceR);
        layerConfig.rminmax.second =
            std::max(layerConfig.rminmax.second, surfaceR);
        // Check for splitting and record the radii
        if (m_cfg.centralLayerSplit > 0.) {
          bool foundR = false;
          for (auto& sradius : layerConfig.splitRadii) {
            if (std::abs(surfaceR - sradius) < m_cfg.centralLayerSplit) {
              foundR = true;
            }
          }
          if (not foundR) {
            layerConfig.splitRadii.push_back(surfaceR);
          }
        }
      }
    } else {
      // is not yet the senstive one
      ACTS_VERBOSE(offset << "[<<] not accepted, stepping down.");
      // set the visibility to kFALSE
      if (m_cfg.setVisibility) {
        tgNode->SetVisibility(kFALSE);
      }
      // screen output for disk debugging
      if (type != 0) {
        ACTS_VERBOSE(offset << "  node translation in z = " << z);
      }
      // build the matrix
      TGeoHMatrix nTransform =
          TGeoCombiTrans(tgTransform) * TGeoCombiTrans(*tgMatrix);
      std::string suffix = "_transform";
      nTransform.SetName((tNodeName + suffix).c_str());
      // if it's not accepted, get the associated volume
      TGeoVolume* nodeVolume = tgNode->GetVolume();
      // step down one further
      resolveSensitive(gctx, layerSurfaces, nodeVolume, nullptr, nTransform,
                       layerConfig, type, correctBranch, offset + "  ");
    }
  } else {
    ACTS_VERBOSE("No node present.");
  }
}
