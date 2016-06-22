// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/TGeoPlugins/TGeoLayerBuilder.hpp"
#include "ACTS/Plugins/TGeoPlugins/TGeoDetectorElement.hpp"
#include "ACTS/Tools/ILayerCreator.hpp"
#include "ACTS/Utilities/MsgMacros.hpp"
#include "TGeoManager.h"
#include <stdio.h>

Acts::TGeoLayerBuilder::TGeoLayerBuilder(
    const Acts::TGeoLayerBuilder::Config& config)
  : m_cfg()
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

const Acts::LayerVector
Acts::TGeoLayerBuilder::negativeLayers() const
{
  LayerVector nVector;
  return std::move(nVector);
}

const Acts::LayerVector
Acts::TGeoLayerBuilder::centralLayers() const
{
  // the return vector
  LayerVector cVector;
  // bail out if you have no gGeoManager
  if (!gGeoManager) return std::move(cVector);

  MSG_INFO("Central Layers : found " << m_cfg.centralLayerConfigs.size()
                                     << " configurations.");

  for (auto layerCfg : m_cfg.centralLayerConfigs) {
    MSG_INFO("- layer configuration found for layer " << layerCfg.layerName
                                                      << " with sensor "
                                                      << layerCfg.sensorName);
    TGeoVolume* volume = gGeoManager->GetVolume(layerCfg.layerName.c_str());
    if (volume) {
      // prepare the vector for the sensitive nodes
      std::vector<TGeoNode*> sensitiveNodes;
      // run recursive collection
      collectSensitive(volume, nullptr, layerCfg.sensorName, sensitiveNodes);
      MSG_INFO("- layer found to have " << sensitiveNodes.size()
                                        << " sensitive sensors ");
      // create the detector surface vector
      std::vector<const Acts::Surface*> detSurfaces;
      detSurfaces.reserve(sensitiveNodes.size());
      // loop and fill
      for (auto& sNode : sensitiveNodes) {
        // create the new detector element
        auto tgElement
            = std::make_shared<Acts::TGeoDetectorElement>(Identifier(), sNode);
        m_elementStore.push_back(tgElement);
        // create a layer out of the surfaces
        detSurfaces.push_back(&(tgElement->surface()));
        // create the layer
        cVector.push_back(m_cfg.layerCreator->cylinderLayer(
            detSurfaces, 1., 5., layerCfg.binsLoc0, layerCfg.binsLoc0));
      }
    }
  }

  return std::move(cVector);
}

const Acts::LayerVector
Acts::TGeoLayerBuilder::positiveLayers() const
{
  LayerVector pVector;
  return std::move(pVector);
}

void
Acts::TGeoLayerBuilder::collectSensitive(
    TGeoVolume*             tgVolume,
    TGeoNode*               tNode,
    const std::string&      sensitiveName,
    std::vector<TGeoNode*>& sensitiveNodes) const
{
  // if it's still the master volume
  if (tgVolume) {
    auto  daugthers = tgVolume->GetNodes();
    TIter iObj(daugthers);
    while (TObject* obj = iObj()) {
      // dynamic_cast to a node
      TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
      if (node) collectSensitive(nullptr, node, sensitiveName, sensitiveNodes);
    }

  } else {
    // get the name as a string and compare
    std::string tNodeName = tNode->GetName();
    MSG_VERBOSE("-- node :" << tNodeName);
    if (tNodeName.find(sensitiveName) != std::string::npos) {
      // senstive volume found, collect it
      sensitiveNodes.push_back(tNode);
    }
    // get the children nodes from the
    auto  daugthers = tNode->GetNodes();
    TIter iObj(daugthers);
    while (TObject* obj = iObj()) {
      // dynamic_cast to a node
      TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
      if (node) collectSensitive(nullptr, node, sensitiveName, sensitiveNodes);
    }
  }
}
