// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialMapping.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Plugins/MaterialPlugins/MaterialMapping.hpp"
#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Plugins/MaterialPlugins/LayerMaterialRecord.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Helpers.hpp"

Acts::MaterialMapping::MaterialMapping(const Config&           cnf,
                                       std::unique_ptr<Logger> log)
  : m_cnf(cnf), m_logger(std::move(log))
{
  // check if extrapolation engine is given
  if (!m_cnf.extrapolationEngine) {
    ACTS_ERROR("[!]No extrapolation engine given!");
  } else
    ACTS_INFO("Extrapolation engine successfully retrieved!");
}

Acts::MaterialMapping::~MaterialMapping()
{
}

void
Acts::MaterialMapping::setLogger(std::unique_ptr<Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

void
Acts::MaterialMapping::mapMaterial(const MaterialTrackRecord& matTrackRec)
{
  // create object which connects layer with hits
  std::vector<std::pair<const Acts::Layer*, Acts::Vector3D>> layersAndHits;
  // associate the m material to the layer only if hits have been collected
  if (collectLayersAndHits(matTrackRec, layersAndHits))
    associateLayerMaterial(matTrackRec, layersAndHits);
}

bool
Acts::MaterialMapping::collectLayersAndHits(
    const MaterialTrackRecord& matTrackRec,
    std::vector<std::pair<const Acts::Layer*, Acts::Vector3D>>& layersAndHits)
{
  // access the parameters
  double                    eta           = matTrackRec.eta();
  double                    phi           = matTrackRec.phi();
  std::vector<MaterialStep> materialSteps = matTrackRec.materialSteps();
  Acts::Vector3D            startPos(matTrackRec.position().x,
                          matTrackRec.position().y,
                          matTrackRec.position().z);
  // let's extrapolate through the detector and remember which layers(with
  // material) should have been hit
  if (materialSteps.size()) {
    // get the number of materialsteps
    ACTS_DEBUG("Successfuly retrieved " << materialSteps.size()
                                        << "materialSteps");
    // propagate through the detector and collect the layers hit in the given
    // direction eta phi
    // calculate the direction in cartesian coordinates
    Acts::Vector3D direction(cos(phi), sin(phi), sinh(eta));
    // create the beginning neutral parameters to extrapolate through the
    // geometry
    std::unique_ptr<Acts::ActsSymMatrixD<5>> cov;
    Acts::NeutralCurvilinearParameters       startParameters(
        std::move(cov), startPos, direction);
    // c        Acts::CurvilinearParameters
    // startParameters(std::move(cov),startPos,direction,1.);
    // create a neutral extrapolation cell and configure it to only collect
    // layer and surfaces with a SurfaceMaterialProxy
    Acts::ExtrapolationCell<Acts::NeutralParameters> ecc(startParameters);
    // c        Acts::ExtrapolationCell<Acts::TrackParameters>
    // ecc(startParameters);
    ecc.addConfigurationMode(Acts::ExtrapolationMode::StopAtBoundary);
    ecc.addConfigurationMode(Acts::ExtrapolationMode::FATRAS);
    ecc.addConfigurationMode(Acts::ExtrapolationMode::CollectSensitive);
    ecc.addConfigurationMode(Acts::ExtrapolationMode::CollectMaterial);

    // call the extrapolation engine
    // screen output
    ACTS_DEBUG("===> forward extrapolation - collecting material layers <<===");
    // call the extrapolation engine
    Acts::ExtrapolationCode eCode = m_cnf.extrapolationEngine->extrapolate(ecc);
    // end parameter, if there
    if (eCode.isSuccess()) {
      // number of layers hit
      size_t nLayersHit = ecc.extrapolationSteps.size();
      ACTS_VERBOSE("[+] Extrapolation to layers did succeed and found "
                   << nLayersHit
                   << " layers.");
      layersAndHits.reserve(nLayersHit);
      // find all the intersected material - remember the last parameters
      std::unique_ptr<const Acts::NeutralParameters> parameters = nullptr;
      // loop over the collected information
      for (auto& es : ecc.extrapolationSteps) {
        if (es.stepConfiguration.checkMode(
                Acts::ExtrapolationMode::CollectMaterial)) {
          // CHECK now different than original, & raw pointer to layer
          layersAndHits.push_back(
              std::make_pair(es.layer, es.materialPosition));
        }
        // continue if we have parameters
      }  // loop over extrapolationsteps
    }    // extrapolation success
  }      // stepCollection
  return bool(layersAndHits.size());
}

void
Acts::MaterialMapping::associateLayerMaterial(
    const MaterialTrackRecord& matTrackRec,
    std::vector<std::pair<const Acts::Layer*, Acts::Vector3D>>& layersAndHits)
{
  // now go through the material step collection and find best fitting layer
  // layers are ordered, hence you can move the starting point along
  // CHECK if Pre-Post-Update needs to be done - Here!
  size_t                    nLayersHit    = layersAndHits.size();
  size_t                    currentLayer  = 0;
  std::vector<MaterialStep> materialSteps = matTrackRec.materialSteps();
  // loop through hits and find the closest layer, the start point moves
  // outwards as we go
  for (auto& step : materialSteps) {
    ACTS_VERBOSE("[L] starting from layer "
                 << currentLayer
                 << " from layer collection for this step.");
    // step length and position
    Acts::Vector3D pos(step.position().x, step.position().y, step.position().z);
    // now find the closest layer
    // if the currentlayer is the last layer and the hit is still inside ->
    // assign & check if the layers before have been assigned the right way -
    // reassign in case another layer fits better
    if (currentLayer < nLayersHit - 1) {
      // search through the layers - this is the reference distance for
      // projection
      double currentDistance
          = (pos - layersAndHits.at(currentLayer).second).mag();
      ACTS_VERBOSE("  - current distance is "
                   << currentDistance
                   << " from "
                   << Acts::toString(pos)
                   << " and "
                   << Acts::toString(layersAndHits.at(currentLayer).second));
      // check if other layer is more suitable
      for (size_t testLayer = (currentLayer + 1); testLayer < nLayersHit;
           ++testLayer) {
        // calculate the distance to the testlayer
        double testDistance = (pos - layersAndHits.at(testLayer).second).mag();
        ACTS_VERBOSE("[L] Testing layer "
                     << testLayer
                     << " from layer collection for this step.");
        ACTS_VERBOSE(" - test distance is "
                     << testDistance
                     << " from "
                     << Acts::toString(pos)
                     << " and "
                     << Acts::toString(layersAndHits.at(testLayer).second));
        if (testDistance < currentDistance) {
          ACTS_VERBOSE(
              "[L] Skipping over to current layer " << testLayer << " because "
                                                    << testDistance
                                                    << " < "
                                                    << currentDistance);
          // the test distance did shrink - update currentlayer
          currentLayer    = testLayer;
          currentDistance = testDistance;
        }  // testdistance < currentDistance
        else
          break;  // stick to the layer you have
      }           // check for better fitting layers
    }             // if last layer
    // the current layer *should* be correct now
    const Acts::Layer* assignedLayer    = layersAndHits.at(currentLayer).first;
    Acts::Vector3D     assignedPosition = layersAndHits.at(currentLayer).second;
    // create material Properties
    const Acts::MaterialProperties* layerMaterialProperties
        = new MaterialProperties(step.material());
    // associate the hit
    ACTS_VERBOSE("[L] Now associate hit at " << Acts::toString(pos));
    associateHit(assignedLayer, assignedPosition, layerMaterialProperties);
  }  // go through material step collection
}

void
Acts::MaterialMapping::associateHit(
    const Layer*                    layer,
    const Acts::Vector3D&           position,
    const Acts::MaterialProperties* layerMaterialProperties)
{
  auto layerRecord = m_layerRecords.find(layer);
  // if layer was not present already create new Material Record
  if (layerRecord == m_layerRecords.end()) {
    // get the bin utility
    const Acts::BinUtility* binUtility = layer->material()->binUtility();
    // create the material record
    ACTS_VERBOSE("[L] Creating new Layer record at position "
                 << Acts::toString(position));
    m_layerRecords[layer] = LayerMaterialRecord(binUtility);
  }
  ACTS_VERBOSE("[L] Add new layer material properties  at position "
               << Acts::toString(position));
  ACTS_VERBOSE(*layerMaterialProperties);
  // add material to record, if record exists already
  m_layerRecords[layer].addLayerMaterialProperties(position,
                                                   layerMaterialProperties);
}

void
Acts::MaterialMapping::averageLayerMaterial()
{
  ACTS_VERBOSE(m_layerRecords.size() << " LayerMaterialRecords to be averaged");
  // average over the layer material
  for (auto& layRecord : m_layerRecords) {
    layRecord.second.averageMaterial();
  }
}

void
Acts::MaterialMapping::finalizeLayerMaterial()
{
  ACTS_VERBOSE(m_layerRecords.size()
               << " LayerMaterialRecords to be finalized");
  // finally set the material of the layers
  for (auto& layRecord : m_layerRecords) {
    layRecord.first->materialSurface()->setAssociatedMaterial(
        layRecord.second.layerMaterial());
  }
}
