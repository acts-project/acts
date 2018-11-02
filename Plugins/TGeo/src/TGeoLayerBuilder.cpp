// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"
#include <stdio.h>
#include "Acts/Layers/ProtoLayer.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

Acts::TGeoLayerBuilder::TGeoLayerBuilder(
    const Acts::TGeoLayerBuilder::Config& config,
    std::unique_ptr<const Logger>         logger)
  : m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(config);
}

Acts::TGeoLayerBuilder::~TGeoLayerBuilder() = default;

void
Acts::TGeoLayerBuilder::setConfiguration(
    const Acts::TGeoLayerBuilder::Config& config)
{
  m_cfg = config;
}

void
Acts::TGeoLayerBuilder::setLogger(std::unique_ptr<const Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

const Acts::LayerVector
Acts::TGeoLayerBuilder::negativeLayers() const
{
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto        mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector nVector;
  mutableThis->buildLayers(nVector, -1);
  return nVector;
}

const Acts::LayerVector
Acts::TGeoLayerBuilder::centralLayers() const
{
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto        mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector cVector;
  mutableThis->buildLayers(cVector, 0);
  return cVector;
}

const Acts::LayerVector
Acts::TGeoLayerBuilder::positiveLayers() const
{
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto        mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector pVector;
  mutableThis->buildLayers(pVector, -1);
  return pVector;
}

void
Acts::TGeoLayerBuilder::buildLayers(LayerVector& layers, int type)
{

  // bail out if you have no gGeoManager
  if (gGeoManager == nullptr) {
    return;
  }

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
    std::vector<std::shared_ptr<const Surface>> layerSurfaces;

    ACTS_DEBUG("- layer configuration found for layer " << layerCfg.layerName
                                                        << " with sensor "
                                                        << layerCfg.sensorName);
    // we have to step down from the top volume each time to collect the logical
    // tree
    TGeoVolume* tvolume = gGeoManager->GetTopVolume();
    if (tvolume != nullptr) {
      // recursively step down
      resolveSensitive(
          layerSurfaces, tvolume, nullptr, TGeoIdentity(), layerCfg, type);
      // screen output
      ACTS_DEBUG(
          "- number of senstive sensors found : " << layerSurfaces.size());
      // create the layer  - either way
      if (type == 0) {
        ProtoLayer pl(layerSurfaces);
        pl.envR = {layerCfg.envelope.first, layerCfg.envelope.second};
        pl.envZ = {layerCfg.envelope.second, layerCfg.envelope.second};
        layers.push_back(m_cfg.layerCreator->cylinderLayer(
            layerSurfaces, layerCfg.binsLoc0, layerCfg.binsLoc1, pl));
      } else {
        ProtoLayer pl(layerSurfaces);
        pl.envR = {layerCfg.envelope.first, layerCfg.envelope.second};
        pl.envZ = {layerCfg.envelope.second, layerCfg.envelope.second};
        layers.push_back(m_cfg.layerCreator->discLayer(
            layerSurfaces, layerCfg.binsLoc0, layerCfg.binsLoc1, pl));
      }
    }
  }
}

void
Acts::TGeoLayerBuilder::resolveSensitive(
    std::vector<std::shared_ptr<const Acts::Surface>>& layerSurfaces,
    TGeoVolume*                                        tgVolume,
    TGeoNode*                                          tgNode,
    const TGeoMatrix&                                  tgTransform,
    const LayerConfig&                                 layerConfig,
    int                                                type,
    bool                                               correctBranch,
    const std::string&                                 offset)
{
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
    if (!correctVolume
        && (volumeName.find(layerConfig.layerName) != std::string::npos
            || match(layerConfig.layerName.c_str(), volumeName.c_str()))) {
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
        resolveSensitive(layerSurfaces,
                         nullptr,
                         node,
                         tgTransform,
                         layerConfig,
                         type,
                         correctVolume,
                         offset + "  ");
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
    bool branchHit
        = correctBranch || (layerConfig.sensorName == layerConfig.layerName);
    if (branchHit
        && (tNodeName.find(layerConfig.sensorName) != std::string::npos
            || match(layerConfig.sensorName.c_str(), tNodeName.c_str()))) {

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
            Identifier(),
            tgNode,
            &tgTransform,
            layerConfig.localAxes,
            m_cfg.unit);
        // record the element @todo solve with provided cache
        m_elementStore.push_back(tgElement);
        // record the surface
        // element owns the surface, we give shared ownership to
        // layer surfaces -> i.e. the Layer to be built
        layerSurfaces.push_back(tgElement->surface().getSharedPtr());
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
      TGeoHMatrix nTransform
          = TGeoCombiTrans(tgTransform) * TGeoCombiTrans(*tgMatrix);
      std::string suffix = "_transform";
      nTransform.SetName((tNodeName + suffix).c_str());
      // if it's not accepted, get the associated volume
      TGeoVolume* nodeVolume = tgNode->GetVolume();
      // step down one further
      resolveSensitive(layerSurfaces,
                       nodeVolume,
                       nullptr,
                       nTransform,
                       layerConfig,
                       type,
                       correctBranch,
                       offset + "  ");
    }
  } else {
    ACTS_VERBOSE("No node present.");
  }
}
