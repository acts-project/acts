// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/TGeoPlugins/TGeoLayerBuilder.hpp"
#include <stdio.h>
#include "ACTS/Plugins/TGeoPlugins/TGeoDetectorElement.hpp"
#include "ACTS/Tools/ILayerCreator.hpp"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

Acts::TGeoLayerBuilder::TGeoLayerBuilder(
    const Acts::TGeoLayerBuilder::Config& config,
    std::unique_ptr<Logger>               logger)
  : m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(config);
}

Acts::TGeoLayerBuilder::~TGeoLayerBuilder()
{
}

void
Acts::TGeoLayerBuilder::setConfiguration(
    const Acts::TGeoLayerBuilder::Config& config)
{
  m_cfg = config;
}

void
Acts::TGeoLayerBuilder::setLogger(std::unique_ptr<Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

const Acts::LayerVector
Acts::TGeoLayerBuilder::negativeLayers() const
{
  LayerVector nVector;
  buildLayers(nVector, -1);
  return std::move(nVector);
}

const Acts::LayerVector
Acts::TGeoLayerBuilder::centralLayers() const
{
  LayerVector cVector;
  buildLayers(cVector, 0);
  return std::move(cVector);
}

const Acts::LayerVector
Acts::TGeoLayerBuilder::positiveLayers() const
{
  LayerVector pVector;
  buildLayers(pVector, -1);
  return std::move(pVector);
}

void
Acts::TGeoLayerBuilder::buildLayers(LayerVector& layers, int type) const
{
  // bail out if you have no gGeoManager
  if (!gGeoManager) return;

  // Prepare which ones to build
  std::vector<LayerConfig> layerConfigs;
  std::string              layerType = "No";
  switch (type) {
  case -1: {
    layerConfigs = m_cfg.negativeLayerConfigs;
    layerType    = "Negative";
  } break;
  case 0: {
    layerConfigs = m_cfg.centralLayerConfigs;
    layerType    = "Central";
  } break;
  case 1: {
    layerConfigs = m_cfg.positiveLayerConfigs;
    layerType    = "Positive";
  } break;
  }
  // screen output
  ACTS_DEBUG(layerType << " Layers : found " << layerConfigs.size()
                       << " configurations.");
  for (auto layerCfg : layerConfigs) {
    // prepare the layer surfaces
    std::vector<const Surface*> layerSurfaces;

    ACTS_DEBUG("- layer configuration found for layer " << layerCfg.layerName
                                                        << " with sensor "
                                                        << layerCfg.sensorName);
    // we have to step down from the top volume each time to collect the logical
    // tree
    TGeoVolume* tvolume = gGeoManager->GetTopVolume();
    if (tvolume) {
      // recursively step down
      collectSensitive(
          layerSurfaces, tvolume, nullptr, TGeoIdentity(), layerCfg, type);
      // screen output
      ACTS_DEBUG(
          "- number of senstive sensors found : " << layerSurfaces.size());
      // create the layer  - either way
      if (!type)
        layers.push_back(
            m_cfg.layerCreator->cylinderLayer(layerSurfaces,
                                              layerCfg.envelope.first,
                                              layerCfg.envelope.second,
                                              layerCfg.binsLoc0,
                                              layerCfg.binsLoc1));
      else
        layers.push_back(m_cfg.layerCreator->discLayer(layerSurfaces,
                                                       layerCfg.envelope.first,
                                                       layerCfg.envelope.first,
                                                       layerCfg.envelope.second,
                                                       layerCfg.binsLoc0,
                                                       layerCfg.binsLoc1));
    }
  }
}

void
Acts::TGeoLayerBuilder::collectSensitive(
    std::vector<const Acts::Surface*>& layerSurfaces,
    TGeoVolume*                        tgVolume,
    TGeoNode*                          tgNode,
    const TGeoMatrix&                  tgTransform,
    const LayerConfig&                 layerConfig,
    int                                type,
    bool                               correctBranch,
    const std::string&                 offset) const
{
  /// some screen output for disk debugging
  if (type) {
    const Double_t* ctranslation = tgTransform.GetTranslation();
    ACTS_VERBOSE(offset << "current z translation is : " << ctranslation[2]);
  }

  if (tgVolume) {
    std::string volumeName = tgVolume->GetName();
    /// some screen output indicating that the volume was found
    ACTS_VERBOSE(offset << "[o] Volume : " << volumeName);
    // once in the current branch, always in the current branch
    bool correctVolume = correctBranch;
    if (correctVolume == false
        && volumeName.find(layerConfig.layerName) != std::string::npos) {
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
      if (node)
        collectSensitive(layerSurfaces,
                         nullptr,
                         node,
                         tgTransform,
                         layerConfig,
                         type,
                         correctVolume,
                         offset + "  ");
    }
  }
  /// if you have a node, get the volume and step down further
  if (tgNode) {
    // get the matrix of the current
    const TGeoMatrix* tgMatrix = tgNode->GetMatrix();
    /// get the translation of the parent
    const Double_t* translation = tgTransform.GetTranslation();
    // get the z value
    double z = translation[2];
    // get the name of the node
    std::string tNodeName = tgNode->GetName();
    ACTS_VERBOSE(offset << "[>] Node : " << tNodeName);
    if (correctBranch
        && tNodeName.find(layerConfig.sensorName) != std::string::npos) {
      // set the visibility to kTrue
      if (m_cfg.setVisibility) tgNode->SetVisibility(kTRUE);
      // create the detector element - check on the type for the size
      if (!type || type * z > 0.) {
        //  senstive volume found, collect it
        ACTS_VERBOSE(offset << "[>>] accepted !");
        // create the element
        auto tgElement
            = std::make_shared<Acts::TGeoDetectorElement>(Identifier(),
                                                          tgNode,
                                                          &tgTransform,
                                                          layerConfig.localAxes,
                                                          m_cfg.unit);
        // record the element @TODO solve with provided cache
        m_elementStore.push_back(tgElement);
        // record the surface
        layerSurfaces.push_back(&(tgElement->surface()));
      }
    } else {
      // is not yet the senstive one
      ACTS_VERBOSE(offset << "[<<] not accepted, stepping down.");
      // set the visibility to kFALSE
      if (m_cfg.setVisibility) tgNode->SetVisibility(kFALSE);
      // screen output for disk debugging
      if (type) ACTS_VERBOSE(offset << "  node translation in z = " << z);
      // build the matrix
      TGeoHMatrix nTransform = (*tgMatrix) * tgTransform;
      std::string suffix     = "_transform";
      nTransform.SetName((tNodeName + suffix).c_str());
      // if it's not accepted, get the associated volume
      TGeoVolume* nodeVolume = tgNode->GetVolume();
      // step down one further
      collectSensitive(layerSurfaces,
                       nodeVolume,
                       nullptr,
                       nTransform,
                       layerConfig,
                       type,
                       correctBranch,
                       offset + "  ");
    }
  }
}
