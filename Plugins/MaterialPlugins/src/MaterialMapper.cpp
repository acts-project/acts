// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialMapper.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Plugins/MaterialPlugins/MaterialMapper.hpp"
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
#include "ACTS/Utilities/GeometryObjectSorter.hpp"
#include "ACTS/Material/SurfaceMaterialProxy.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Material/BinnedSurfaceMaterial.hpp"
#include <climits>

Acts::MaterialMapper::MaterialMapper(const Config&           cfg,
                                       std::unique_ptr<Logger> log)
  : m_cfg(cfg), m_logger(std::move(log))
{
  // check if extrapolation engine is given
  if (!m_cfg.extrapolationEngine) {
    ACTS_ERROR("[!]No extrapolation engine given!");
  } else
    ACTS_INFO("Extrapolation engine successfully retrieved!");
}

Acts::MaterialMapper::~MaterialMapper()
{
}

void
Acts::MaterialMapper::setLogger(std::unique_ptr<Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

Acts::MaterialMapper::Cache
Acts::MaterialMapper::materialMappingCache(
    const TrackingGeometry& tGeometry) const
{
  // parse the geometry and find all surfaces with material proxies
  auto world = tGeometry.highestTrackingVolume();
  // create the map
  std::map<GeometryID, SurfaceMaterialRecord> sMap;
  // fill it
  collectMaterialSurfaces(sMap, *world);
  
  ACTS_DEBUG(sMap.size() <<  " Surfaces with PROXIES collected ... ");
  for (auto& smg : sMap) {
    //print out the information
    size_t volumeID  = smg.first.value(GeometryID::volume_mask);
    size_t layerID   = smg.first.value(GeometryID::layer_mask);
    ACTS_VERBOSE(" -> Surface in volume " << volumeID 
                         << " for layer " << layerID);
  }  
  // return it
  return MaterialMapper::Cache(sMap);
}

bool
Acts::MaterialMapper::mapMaterialTrackRecord(
    Cache&                     mappingCache,
    const MaterialTrackRecord& trackRecord) const
{
  // access the parameters
  double   theta = trackRecord.theta();
  double   phi   = trackRecord.phi();
  auto     spos  = trackRecord.position();
  Vector3D vertex(spos.x, spos.y, spos.z);

  // counter for track reconrds
  mappingCache.trackRecordCounter++;

  // get the steps from the detailed geometry
  std::vector<MaterialStep> materialSteps = trackRecord.materialSteps();

  // get the steps from the tracking geometry
  std::vector<AssignedMaterialSteps> assignedSteps;

  // let's extrapolate through the ACTS detector and record all surfaces
  // that have a material proxy
  if (materialSteps.size()) {
    // get the number of materialsteps
    ACTS_DEBUG("Successfuly retrieved " << materialSteps.size()
                                        << " materialSteps");
    // propagate through the detector and collect the layers hit in the given
    // direction eta phi
    // calculate the direction in cartesian coordinates
    Vector3D direction(
        cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    // create the beginning neutral parameters to extrapolate through the
    // geometry
    NeutralCurvilinearParameters startParameters(nullptr, vertex, direction);
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
    ACTS_DEBUG("===> forward extrapolation - collecting intersections <<===");
    // call the extrapolation engine
    ExtrapolationCode eCode = m_cfg.extrapolationEngine->extrapolate(ecc);
    // end parameter, if there
    if (eCode.isSuccess()) {
      // number of layers hit
      size_t nSurfacesHit = ecc.extrapolationSteps.size();
      ACTS_VERBOSE("[+] Extrapolation succeeded resulting in "
                   << nSurfacesHit
                   << " intersected surfaces.");
      // loop over the collected information
      for (auto& es : ecc.extrapolationSteps) {
        if (es.configuration.checkMode(ExtrapolationMode::CollectMaterial)) {
          // geo ID from the surface
          GeometryID assignID = es.surface->geoID();
          // more screen output for debugging
          ACTS_VERBOSE("Material proxy found on surface with ID " << assignID.value());
          // collect the assigned ones
          assignedSteps.push_back(AssignedMaterialSteps(assignID, es.position));
        }
      } // loop over extrapolationsteps
    } // extrapolation success

    // now check how many have been assigned
    ACTS_VERBOSE("[+] Selected " << assignedSteps.size() << " for mapping.")
    
  } // stepCollection

  // run the step assignment
  assignSteps(materialSteps, assignedSteps);

  // and now we fill it into the record
  for (auto& aSteps : assignedSteps) {
    /// panic if the assignedGeoID is not in the mappingCache
    if (mappingCache.surfaceMaterialRecords.find(aSteps.assignedGeoID)
        == mappingCache.surfaceMaterialRecords.end()){
      // that deserves a WARNING
      ACTS_WARNING("[-] Material surface with " 
                   <<  aSteps.assignedGeoID.value() << " not found in cache.");
      continue;
      
    } else
    mappingCache.surfaceMaterialRecords[aSteps.assignedGeoID].assignMaterialSteps(aSteps);
  }
  
  return true;
}

std::map<Acts::GeometryID, Acts::SurfaceMaterial*> 
Acts::MaterialMapper::createSurfaceMaterial(Cache& mappingCache) const 
{
  // the return map for the surface material
  std::map<GeometryID, SurfaceMaterial*> surfaceMaterialMap;
  ACTS_DEBUG("Creating material maps for surfaces");
  // let's loop over the surface material records and create the corresponding maps
  for (auto& smr : mappingCache.surfaceMaterialRecords){
    // get the corresponding GeometryID
    GeometryID surfaceID = smr.first;
    // some more screen output
    ACTS_DEBUG(" -> Volume | Layer | Approach | Sensitive  "
               << surfaceID.value(GeometryID::volume_mask) << " | "
               << surfaceID.value(GeometryID::layer_mask)  << " | "
               << surfaceID.value(GeometryID::approach_mask) << " | "
               << surfaceID.value(GeometryID::sensitive_mask));
    // get the BinUtility
    const BinUtility& bUtility = smr.second.binUtility();
    // get the mapped material
    const MaterialRecord& mMaterial = smr.second.mappedMaterial();
    // the size of the map to be created
    size_t bins0 = bUtility.bins(0);
    size_t bins1 = bUtility.bins(1);
    // even more screen output
    ACTS_VERBOSE(" -> Material matrix with [ " << bins0 << " x " <<  bins1 << " ] bins" );
    // prepare the matrix
    MaterialPropertiesVector mVector(bins0,nullptr);
    MaterialPropertiesMatrix mMatrix(bins1,mVector);
    // fill the matrix
    for (size_t i1 = 0; i1 < bins1; ++i1){
      for (size_t i0 = 0; i0 < bins0; ++i0 ){
        // default is nullptr
        MaterialProperties* binMaterial = nullptr;
        // get what you have
        size_t mappingHits = mMaterial[i1][i0].second;
        if (mappingHits){
          // the statistical scalor
          double statScalor = 1./double(mappingHits);
          // very verbose screen output
          ACTS_VERBOSE(" --[ " << i0 << " x " << i1 << " ] has " << mappingHits << " entries");
          // take the mapped material 
          MaterialProperties matProperties = mMaterial[i1][i0].first;
          double thickness = statScalor*matProperties.thickness();
          double rho       = statScalor*matProperties.averageRho();
          double A         = statScalor*matProperties.averageA();
          double Z         = statScalor*matProperties.averageZ();
          double tInX0     = statScalor*matProperties.thicknessInX0();
          double tInL0     = statScalor*matProperties.thicknessInL0();
          // recreate X0, L0 
          float x0 = (thickness != 0. && tInX0 != 0.) ? thickness / tInX0 : 0.;
          float l0 = (thickness != 0. && tInL0 != 0.) ? thickness / tInL0 : 0.;
          // create the new material (with statistical weights)
          binMaterial = new MaterialProperties(x0, l0, A, Z, rho, thickness);
        }
        mMatrix[i1][i0] = binMaterial;
      }
    }
    // create surface material
    surfaceMaterialMap[surfaceID] 
      = new BinnedSurfaceMaterial(bUtility,
                                  mMatrix,
                                  0.,
                                  mappingCache.trackRecordCounter);
  }
  // now return the material map  
  return surfaceMaterialMap;
}

void
Acts::MaterialMapper::assignSteps(
    const std::vector<MaterialStep>& materialSteps,
    std::vector<AssignedMaterialSteps>& assignedSteps) const
{

  // we will rely on the fact that the material steps are ordered
  // and so are the assigned steps

  // the iterators
  std::vector<AssignedMaterialSteps>::iterator asIter      = assignedSteps.begin();
  std::vector<AssignedMaterialSteps>::iterator asIterFlush = assignedSteps.begin();
  std::vector<AssignedMaterialSteps>::iterator asIterLast  = assignedSteps.end()-1;
  std::vector<AssignedMaterialSteps>::iterator asIterEnd   = assignedSteps.end();

  // loop over the steps
  for (auto& mStep : materialSteps) {
    // start with infinity distance
    double lDist2 = std::numeric_limits<double>::infinity();
    // the material position as a Vector3D
    Vector3D mPosition(
        mStep.position().x, mStep.position().y, mStep.position().z);
    // now assign to a step
    asIter = asIterFlush;
    for (; asIter != asIterEnd; ++asIter) {
      // the currentDist
      double cDist2 = (mPosition - asIter->assignedPosition).mag2();
      if (cDist2 < lDist2) {
        // set the new closest distance
        lDist2 = cDist2;
        // remember where you are
        asIterFlush = asIter;
        if (asIterFlush == asIterLast){
          // the last iteration point wins the step
          asIterFlush->assignedSteps.push_back(mStep);
          // just to avoid the next if
          break;
        }
      } else if (asIter != assignedSteps.begin()) {
        // distance increase detected and it is not the first step
        // the last iteration point wins the step
        asIterFlush->assignedSteps.push_back(mStep);
        // we set the new start iterator to be the Flush iterator
        asIter = asIterFlush;
        // and break the inner loop, this jumps to the next material
        // step which will try assignment from the new start position
        break;
      }
    }
  }
}

void
Acts::MaterialMapper::collectMaterialSurfaces(
    std::map<GeometryID, SurfaceMaterialRecord>& sMap,
    const TrackingVolume& tVolume) const
{
  
  ACTS_VERBOSE("Checking volume '" << tVolume.volumeName() 
               << "' for material surfaces.")
  // check the boundary surfaces
  for (auto& bSurface : tVolume.boundarySurfaces())
    checkAndInsert(sMap, bSurface->surfaceRepresentation());

  // check the confined layers
  if (tVolume.confinedLayers()) {
    for (auto& cLayer : tVolume.confinedLayers()->arrayObjects()) {
      // take only layers that are not navigation layers
      if (cLayer->layerType() != navigation) {
        // check the representing surface
        checkAndInsert(sMap, cLayer->surfaceRepresentation());
        // get the approach surfaces if present
        if (cLayer->approachDescriptor()) {
          for (auto& aSurface :
               cLayer->approachDescriptor()->containedSurfaces())
            if (aSurface) checkAndInsert(sMap, *aSurface);
        }
        // get the sensitive surface is present
        if (cLayer->surfaceArray()) {
          // sensitive surface loop
          for (auto& sSurface : cLayer->surfaceArray()->arrayObjects())
            if (sSurface) checkAndInsert(sMap, *sSurface);
        }
      }
    }
  }

  // step down into the sub volume
  if (tVolume.confinedVolumes()) {
    for (auto& sVolume : tVolume.confinedVolumes()->arrayObjects()) {
      // recursive call
      collectMaterialSurfaces(sMap, *sVolume);
    }
  }
}

void
Acts::MaterialMapper::checkAndInsert(
    std::map<GeometryID, SurfaceMaterialRecord>& sMap,
    const Surface& surface) const
{

  // check if the surface has a proxy
  if (surface.associatedMaterial()) {
    // we need a dynamic_cast to a surface material proxy
    auto smp = dynamic_cast<const SurfaceMaterialProxy*>(
        surface.associatedMaterial());
    if (smp){
      // get the geo id
      auto geoID = surface.geoID();
      size_t volumeID = geoID.value(GeometryID::volume_mask);
      ACTS_VERBOSE("Material surface found with volumeID " << volumeID);
      ACTS_VERBOSE("       - surfaceID is " << geoID.value());
      sMap[geoID]
          = SurfaceMaterialRecord(surface, *smp->binUtility());

      }
  }
}

// void
// Acts::MaterialMapper::associateLayerMaterial(
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
//  // create object which connects layer with the original material step and
//  its
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
//    Acts::Vector3D pos(step.position().x, step.position().y,
//    step.position().z);
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
//      for (std::vector<std::pair<const Acts::Layer*,
//      Acts::Vector3D>>::iterator
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
//                                 step.material().thickness() /
//                                 pathCorrection);
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
/// Acts::MaterialMapper::associateHit(
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
/// Acts::MaterialMapper::averageLayerMaterial()
/// {
///   ACTS_VERBOSE(m_surfaceMaterialRecords.size() << " SurfaceMaterialRecords
///   to be averaged");
///   // average over the layer material
///   for (auto& layRecord : m_surfaceMaterialRecords) {
///     layRecord.second.averageMaterial();
///   }
/// }
///
/// void
/// Acts::MaterialMapper::finalizeLayerMaterial()
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
