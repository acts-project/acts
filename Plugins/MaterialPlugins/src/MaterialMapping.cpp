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
#include "ACTS/Plugins/MaterialPlugins/SurfaceMaterialRecord.hpp"
#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Helpers.hpp"
#include "ACTS/Material/SurfaceMaterialProxy.hpp"

Acts::MaterialMapping::MaterialMapping(const Config&           cfg,
                                       std::unique_ptr<Logger> log)
  : m_cfg(cfg)
  , m_logger(std::move(log))
{
  // check if extrapolation engine is given
  if (!m_cfg.extrapolationEngine) {
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

std::map<Acts::GeometryID, Acts::SurfaceMaterialRecord >
Acts::MaterialMapping::materialMappingCache(const TrackingGeometry& tGeometry) const
{
  // parse the geometry and find all surfaces with material proxies
  auto world = tGeometry.highestTrackingVolume();
  // create the map
  std::map<Acts::GeometryID, SurfaceMaterialRecord > sMap;
  // fill it
  collectMaterialSurfaces(sMap,*world);
  // return it
  return sMap;
}


bool
Acts::MaterialMapping::mapMaterialTrackRecord(Cache& mappingCache,
                                              const MaterialTrackRecord& trackRecord) const
{
  // access the parameters
  double   theta  = trackRecord.theta();
  double   phi    = trackRecord.phi();
  auto     spos   = trackRecord.position();
  Vector3D vertex(spos.x,spos.y,spos.z);
  
  // get the steps from the detailed geometry
  std::vector<MaterialStep> materialSteps = trackRecord.materialSteps();
  
  // get the steps from the tracking geometry
  std::vector< std::pair< GeometryID, Vector3D > > assignSteps;
  
  // let's extrapolate through the ACTS detector and record all surfaces
  // that have a material proxy
  if (materialSteps.size()) {
    // get the number of materialsteps
    ACTS_DEBUG("Successfuly retrieved " << materialSteps.size()
                                        << "materialSteps");
    // propagate through the detector and collect the layers hit in the given
    // direction eta phi
    // calculate the direction in cartesian coordinates
    Vector3D direction(cos(phi)*sin(theta), sin(phi)*sin(theta),cos(theta));
    // create the beginning neutral parameters to extrapolate through the
    // geometry
    // std::unique_ptr<ActsSymMatrixD<NGlobalPars>> cov;
    Acts::NeutralCurvilinearParameters startParameters(
        nullptr, vertex, direction);
    // create a neutral extrapolation cell and configure it:
    // - to collect surfaces with a SurfaceMaterialProxy
    // - to step at the detector boundary 
    // - to run in a FATRAS type approach 
    ExtrapolationCell<NeutralParameters> ecc(startParameters);
    ecc.addConfigurationMode(ExtrapolationMode::StopAtBoundary);
    ecc.addConfigurationMode(ExtrapolationMode::FATRAS);
    ecc.addConfigurationMode(ExtrapolationMode::CollectSensitive);
    ecc.addConfigurationMode(ExtrapolationMode::CollectMaterial);
    // call the extrapolation engine
    // screen output
    ACTS_DEBUG("===> forward extrapolation - collecting material layers <<===");
    // call the extrapolation engine
    ExtrapolationCode eCode = m_cfg.extrapolationEngine->extrapolate(ecc);
    // end parameter, if there
    if (eCode.isSuccess()) {
      // number of layers hit
      size_t nSurfacesHit = ecc.extrapolationSteps.size();
      ACTS_VERBOSE("[+] Extrapolation to layers did succeed and found "
                   << nSurfacesHit
                   << " layers.");
      // loop over the collected information
      for (auto& es : ecc.extrapolationSteps) {
        if (es.configuration.checkMode(ExtrapolationMode::CollectMaterial)){
          GeometryID      assignID = es.parameters->referenceSurface().geoID();
          const Vector3D& position = es.parameters->position();
          // collect the assigned ones
          assignSteps.push_back(std::pair<GeometryID,Vector3D>(assignID,position));
        }
        // continue if we have parameters
      }  // loop over extrapolationsteps
    }    // extrapolation success
  }      // stepCollection
  return true;
}

void
Acts::MaterialMapping::collectMaterialSurfaces(std::map<GeometryID, SurfaceMaterialRecord>& sMap,
                                               const TrackingVolume& tVolume) const
{
  // check the boundary surfaces
  for (auto& bSurface : tVolume.boundarySurfaces() )
    checkAndInsert(sMap, bSurface->surfaceRepresentation());
  
  // check the confined layers
  if ( tVolume.confinedLayers() ){
    for (auto& cLayer : tVolume.confinedLayers()->arrayObjects()){
      // take only layers that are not navigation layers
      if (cLayer->layerType() != navigation){
        // check the representing surface
        checkAndInsert(sMap, cLayer->surfaceRepresentation());
        // get the approach surfaces if present
        if (cLayer->approachDescriptor()){
          for (auto& aSurface : cLayer->approachDescriptor()->containedSurfaces())
            if (aSurface)
              checkAndInsert(sMap,*aSurface);
        }
        // get the sensitive surface is present
        if (cLayer->surfaceArray()){
          // sensitive surface loop
          for (auto& sSurface : cLayer->surfaceArray()->arrayObjects())
            if(sSurface)
              checkAndInsert(sMap, *sSurface);
        }
      }
    }
  }
  
  // step down into the sub volume
  if (tVolume.confinedVolumes()){
    for (auto& sVolume : tVolume.confinedVolumes()->arrayObjects() ){
      // recursive call
      collectMaterialSurfaces(sMap,*sVolume);
    }
  }
}

void
Acts::MaterialMapping::checkAndInsert(std::map<GeometryID, SurfaceMaterialRecord>& sMap,
                                      const Surface& surface) const
{
  
  // check if the surface has a proxy
  if (surface.associatedMaterial()){
    // we need a dynamic_cast to a surface material proxy
    const SurfaceMaterialProxy* smp
    = dynamic_cast<const SurfaceMaterialProxy*>(surface.associatedMaterial());
    if (smp)
      sMap[surface.geoID()] = SurfaceMaterialRecord(surface, *smp->binUtility());
  }
}



//void
//Acts::MaterialMapping::associateLayerMaterial(
//    const MaterialTrackRecord& trackRecord,
//    std::vector<std::pair<const Acts::Layer*, Acts::Vector3D>>& layersAndHits)
//{
//  // now go through the material step collection and find best fitting layer
//  // layers are ordered, hence you can move the starting point along
//  std::vector<std::pair<const Acts::Layer*, Acts::Vector3D>>::iterator
//      currentLayer
//      = layersAndHits.begin();
//  // access the material steps of this track record
//  std::vector<MaterialStep> materialSteps = trackRecord.materialSteps();
//
//  // create object which connects layer with the original material step and its
//  // assigned position on the layer
//  std::map<const Acts::Layer*,
//           std::pair<const Vector3D, std::vector<MaterialStep>>>
//      layersPosAndSteps;
//
//  // loop through hits and find the closest layer, the start point moves
//  // outwards as we go
//  for (auto& step : materialSteps) {
//    ACTS_VERBOSE("[L] starting from layer "
//                 << std::distance(layersAndHits.begin(), currentLayer)
//                 << " from layer collection for this step.");
//    // step length and position
//    Acts::Vector3D pos(step.position().x, step.position().y, step.position().z);
//    // now find the closest layer
//    // if the currentlayer is the last layer and the hit is still inside ->
//    // assign & check if the layers before have been assigned the right way -
//    // reassign in case another layer fits better
//    if (currentLayer != std::prev(layersAndHits.end())) {
//      // search through the layers - this is the reference distance for
//      // projection
//      double currentDistance = (pos - currentLayer->second).mag();
//      ACTS_VERBOSE(
//          "  - current distance is " << currentDistance << " from "
//                                     << Acts::toString(pos)
//                                     << " and "
//                                     << Acts::toString(currentLayer->second));
//      // check if other layer is more suitable
//      for (std::vector<std::pair<const Acts::Layer*, Acts::Vector3D>>::iterator
//               testLayer
//           = std::next(currentLayer);
//           testLayer != layersAndHits.end();
//           ++testLayer) {
//        // calculate the distance to the testlayer
//        double testDistance = (pos - testLayer->second).mag();
//        ACTS_VERBOSE("[L] Testing layer "
//                     << std::distance(layersAndHits.begin(), testLayer)
//                     << " from layer collection for this step.");
//        ACTS_VERBOSE(
//            " - test distance is " << testDistance << " from "
//                                   << Acts::toString(pos)
//                                   << " and "
//                                   << Acts::toString(testLayer->second));
//        if (testDistance < currentDistance) {
//          ACTS_VERBOSE("[L] Skipping over to current layer "
//                       << std::distance(layersAndHits.begin(), testLayer)
//                       << " because "
//                       << testDistance
//                       << " < "
//                       << currentDistance);
//          // the test distance did shrink - update currentlayer
//          currentLayer    = testLayer;
//          currentDistance = testDistance;
//        }  // testdistance < currentDistance
//        else
//          break;  // stick to the layer you have
//      }           // check for better fitting layers
//    }             // if last layer
//    // the current layer *should* be correct now
//    const Acts::Layer* assignedLayer = currentLayer->first;
//    // correct material thickness with pathcorrection
//    double theta = trackRecord.theta();
//    double phi   = trackRecord.phi();
//    // calculate the direction in cartesian coordinates
//    Acts::Vector3D direction(
//        cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
//    // access the path correction of the associated material surface
//    double pathCorrection
//        = assignedLayer->materialSurface()->pathCorrection(pos, direction);
//
//    // create material Properties
//    const Acts::MaterialProperties* layerMaterialProperties
//        = new MaterialProperties(step.material().material(),
//                                 step.material().thickness() / pathCorrection);
//    // correct also the thickness of the material step
//    Acts::MaterialStep updatedStep(*layerMaterialProperties, step.position());
//    // fill the current material step and its assigned position
//    // first check if layer is already there
//    auto layerPosSteps = layersPosAndSteps.find(assignedLayer);
//    // just fill in material step if layer is already there
//    // otherwise create new entry
//    if (layerPosSteps != layersPosAndSteps.end())
//      layerPosSteps->second.second.push_back(updatedStep);
//    else
//      layersPosAndSteps.emplace(
//          assignedLayer,
//          std::make_pair(currentLayer->second,
//                         std::vector<MaterialStep>{updatedStep}));
//
//    // associate the hit
//    ACTS_VERBOSE("[L] Now associate hit " << Acts::toString(pos) << " at "
//                                          << currentLayer->second);
//  }
//
//  // associate the steps
//  for (auto& layer : layersPosAndSteps) {
//    associateHit(layer.first, layer.second.first, layer.second.second);
//  }
//}

/// void
/// Acts::MaterialMapping::associateHit(
///     const Layer*                           layer,
///     const Acts::Vector3D&                  position,
///     const std::vector<Acts::MaterialStep>& layerMaterialSteps)
/// {
///   auto layerRecord = m_surfaceMaterialRecords.find(layer);
///   // if layer was not present already create new Material Record
///   if (layerRecord == m_surfaceMaterialRecords.end()) {
///     // get the bin utility
///     const Acts::BinUtility* binUtility = layer->material()->binUtility();
///     // create the material record
///     ACTS_VERBOSE("[L] Creating new Layer record, for layer  "
///                  << layer->geoID()
///                  << " at position "
///                  << Acts::toString(position));
///     m_surfaceMaterialRecords[layer] = SurfaceMaterialRecord(binUtility);
///   }
///   ACTS_VERBOSE("[L] Add new layer material properties  at position "
///                << Acts::toString(position));
///   // add material to record, if record exists already
///   m_surfaceMaterialRecords[layer].addLayerMaterialProperties(position,
///                                                    layerMaterialSteps);
/// }
/// 
/// void
/// Acts::MaterialMapping::averageLayerMaterial()
/// {
///   ACTS_VERBOSE(m_surfaceMaterialRecords.size() << " SurfaceMaterialRecords to be averaged");
///   // average over the layer material
///   for (auto& layRecord : m_surfaceMaterialRecords) {
///     layRecord.second.averageMaterial();
///   }
/// }
/// 
/// void
/// Acts::MaterialMapping::finalizeLayerMaterial()
/// {
///   ACTS_VERBOSE(m_surfaceMaterialRecords.size()
///                << " SurfaceMaterialRecords to be finalized");
///   // finally set the material of the layers
///   for (auto& layRecord : m_surfaceMaterialRecords) {
///     // @todo check with Julia how to fix this
///     //layRecord.first->materialSurface()->setAssociatedMaterial(
///     //    layRecord.second.layerMaterial());
///   }
/// }
