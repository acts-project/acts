// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialMapper.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Plugins/MaterialPlugins/MaterialMapper.hpp"
#include <climits>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolation/ExtrapolationCell.hpp"
#include "Acts/Layers/Layer.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Material/SurfaceMaterialProxy.hpp"
#include "Acts/Plugins/MaterialPlugins/SurfaceMaterialRecord.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/GeometryObjectSorter.hpp"
#include "Acts/Utilities/Helpers.hpp"

Acts::MaterialMapper::MaterialMapper(const Config&                 cfg,
                                     std::unique_ptr<const Logger> log)
  : m_cfg(cfg), m_logger(std::move(log))
{
  // check if extrapolation engine is given
  if (!m_cfg.extrapolationEngine) {
    ACTS_ERROR("[!]No extrapolation engine given!");
  } else
    ACTS_DEBUG("Extrapolation engine successfully retrieved!");
}

Acts::MaterialMapper::~MaterialMapper()
{
}

void
Acts::MaterialMapper::setLogger(std::unique_ptr<const Logger> newLogger)
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

  ACTS_DEBUG(sMap.size() << " Surfaces with PROXIES collected ... ");
  for (auto& smg : sMap) {
    // print out the information
    size_t volumeID = smg.first.value(GeometryID::volume_mask);
    size_t layerID  = smg.first.value(GeometryID::layer_mask);
    ACTS_VERBOSE(" -> Surface in volume " << volumeID << " for layer "
                                          << layerID);
  }
  // return it
  return MaterialMapper::Cache(sMap);
}

Acts::MaterialTrack
Acts::MaterialMapper::mapMaterialTrack(Cache&               mappingCache,
                                       const MaterialTrack& materialTrack) const
{
  // access the parameters
  double theta = materialTrack.theta();
  double phi   = materialTrack.phi();
  auto   spos  = materialTrack.position();

  // counter for track reconrds
  mappingCache.materialTrackCounter++;

  // get the steps from the detailed geometry
  std::vector<MaterialStep> materialSteps = materialTrack.materialSteps();

  // get the steps from the tracking geometry
  std::vector<AssignedMaterialSteps> assignedSteps;

  // collect for validation output
  std::map<geo_id_value, std::pair<double, double>> collectedMaterial;

  // let's extrapolate through the Acts detector and record all surfaces
  // that have a material proxy
  if (materialSteps.size()) {
    // get the number of materialsteps
    ACTS_VERBOSE(">> Retrieved " << materialSteps.size()
                                 << " materialSteps along the track at input.");
    // propagate through the detector and collect the layers hit in the given
    // direction eta phi
    // calculate the direction in cartesian coordinates
    Vector3D direction(
        cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
    // create the beginning neutral parameters to extrapolate through the
    // geometry
    NeutralCurvilinearParameters startParameters(nullptr, spos, direction);
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
    // call the extrapolation engine
    ExtrapolationCode eCode = m_cfg.extrapolationEngine->extrapolate(ecc);
    // end parameter, if there
    if (eCode.isSuccess()) {
      // loop over the collected information
      for (auto& es : ecc.extrapolationSteps) {
        if (es.configuration.checkMode(ExtrapolationMode::CollectMaterial)) {
          // geo ID from the surface
          GeometryID assignID = es.surface->geoID();
          // collect the assigned ones
          assignedSteps.push_back(AssignedMaterialSteps(assignID, es.position));
        }
      }  // loop over extrapolationsteps
    }    // extrapolation success

    // now check how many have been assigned
    ACTS_VERBOSE("<< to be mapped onto " << assignedSteps.size()
                                         << " surfaces along the track.");
  }  // stepCollection

  // now associate material steps and assigned steps
  assignSteps(materialSteps, assignedSteps);

  // prepare the return values
  double                    tottX0 = 0.;
  double                    tottL0 = 0.;
  std::vector<MaterialStep> mappedSteps;
  geo_id_value              geoID = 0;

  // and now we fill it into the record
  for (auto& aSteps : assignedSteps) {
    /// panic if the assignedGeoID is not in the mappingCache
    if (mappingCache.surfaceMaterialRecords.find(aSteps.assignedGeoID)
        == mappingCache.surfaceMaterialRecords.end()) {
      // that deserves a WARNING
      ACTS_WARNING("[-] Material surface with " << aSteps.assignedGeoID.value()
                                                << " not found in cache.");
      continue;

    } else {

      // get the surface material record
      auto& surfaceMaterialRecord
          = mappingCache.surfaceMaterialRecords[aSteps.assignedGeoID];
      // the surface to be associated
      auto& surface = surfaceMaterialRecord.surface();
      // get the geo id
      geoID = aSteps.assignedGeoID.value();

      // assign the material
      Vector3D direction(aSteps.assignedPosition.unit());
      double pathC = surface.pathCorrection(aSteps.assignedPosition, direction);
      // loop over the steps and average the material
      // and collpas into one
      double tX0       = 0.;
      double tL0       = 0.;
      double A         = 0.;
      double Z         = 0.;
      double rho       = 0.;
      double thickness = 0.;
      for (auto& currentStep : aSteps.assignedSteps) {
        // thickness and density
        float t       = currentStep.materialProperties().thickness();
        float density = currentStep.materialProperties().averageRho();
        // sum it up
        thickness += t;
        rho += density * t;
        A += currentStep.materialProperties().averageA() * density * t;
        Z += currentStep.materialProperties().averageZ() * density * t;
        // add the thickness in x0
        tX0 += currentStep.materialProperties().thicknessInX0();
        tL0 += currentStep.materialProperties().thicknessInL0();
      }
      // add up the material
      tottX0 += tX0;
      tottL0 += tL0;
      // checks before normalisation
      A /= (rho != 0.) ? rho : 1.;
      Z /= (rho != 0.) ? rho : 1.;
      rho /= (thickness != 0.) ? thickness : 1.;

      // x0 and l0
      double x0 = (thickness != 0. && tX0 != 0.) ? thickness / tX0 : 0.;
      double l0 = (thickness != 0. && tL0 != 0.) ? thickness / tL0 : 0.;

      // create
      // (a) the material and
      Material cMaterial(x0, l0, A, Z, rho);
      // (b) the material properties
      MaterialProperties cMaterialProperties(cMaterial, thickness);
      // (c) the material step
      MaterialStep cMaterialStep(
          cMaterialProperties, aSteps.assignedPosition, geoID);
      // assign the steps to the surface material reccord
      surfaceMaterialRecord.assignMaterialStep(cMaterialStep, pathC);
      // push back the material step
      mappedSteps.push_back(cMaterialStep);
    }
  }
  // return the mapped view of the MaterialTrack
  return MaterialTrack(
      materialTrack.position(), theta, phi, mappedSteps, tottX0, tottL0);
}

std::map<Acts::GeometryID, Acts::SurfaceMaterial*>
Acts::MaterialMapper::createSurfaceMaterial(Cache& mappingCache) const
{
  // the return map for the surface material
  std::map<GeometryID, SurfaceMaterial*> surfaceMaterialMap;
  ACTS_DEBUG("Creating material maps for surfaces");
  // let's loop over the surface material records and create the corresponding
  // maps
  for (auto& smr : mappingCache.surfaceMaterialRecords) {
    // get the corresponding GeometryID
    GeometryID surfaceID = smr.first;
    // some more screen output
    ACTS_DEBUG(" -> Volume | Layer | Approach | Sensitive  "
               << surfaceID.value(GeometryID::volume_mask)
               << " | "
               << surfaceID.value(GeometryID::layer_mask)
               << " | "
               << surfaceID.value(GeometryID::approach_mask)
               << " | "
               << surfaceID.value(GeometryID::sensitive_mask));
    // get the BinUtility
    const BinUtility& bUtility = smr.second.binUtility();
    // get the mapped material
    const MaterialRecord& mMaterial = smr.second.mappedMaterial();
    // the size of the map to be created
    size_t bins0 = bUtility.bins(0);
    size_t bins1 = bUtility.bins(1);
    // even more screen output
    ACTS_VERBOSE(" -> Material matrix with [ " << bins0 << " x " << bins1
                                               << " ] bins");
    // prepare the matrix
    MaterialPropertiesVector mVector(bins0, nullptr);
    MaterialPropertiesMatrix mMatrix(bins1, mVector);
    // fill the matrix
    for (size_t i1 = 0; i1 < bins1; ++i1) {
      for (size_t i0 = 0; i0 < bins0; ++i0) {
        // default is nullptr
        MaterialProperties* binMaterial = nullptr;
        // get what you have
        size_t mappingHits = mMaterial[i1][i0].second;
        if (mappingHits) {
          // the statistical scalor
          double statScalor = 1. / double(mappingHits);
          // very verbose screen output
          ACTS_VERBOSE(" --[ " << i0 << " x " << i1 << " ] has " << mappingHits
                               << " entries");
          // take the mapped material
          MaterialProperties matProperties = mMaterial[i1][i0].first;
          double             thickness = statScalor * matProperties.thickness();
          double             rho = statScalor * matProperties.averageRho();
          double             A   = statScalor * matProperties.averageA();
          double             Z   = statScalor * matProperties.averageZ();
          double             x0  = statScalor * matProperties.material().X0();
          double             l0  = statScalor * matProperties.material().L0();
          // create the new material (with statistical weights)
          binMaterial = new MaterialProperties(x0, l0, A, Z, rho, thickness);
        }
        mMatrix[i1][i0] = binMaterial;
      }
    }
    // create surface material
    surfaceMaterialMap[surfaceID] = new BinnedSurfaceMaterial(
        bUtility, mMatrix, 0., mappingCache.materialTrackCounter);
  }
  // now return the material map
  return surfaceMaterialMap;
}

void
Acts::MaterialMapper::assignSteps(
    const std::vector<MaterialStep>&    materialSteps,
    std::vector<AssignedMaterialSteps>& assignedSteps) const
{

  // will rely on the fact that the
  // material steps are ordered
  // and so are the assigned steps

  // the iterators
  std::vector<AssignedMaterialSteps>::iterator asIter = assignedSteps.begin();
  std::vector<AssignedMaterialSteps>::iterator asIterFlush
      = assignedSteps.begin();
  std::vector<AssignedMaterialSteps>::iterator asIterLast
      = assignedSteps.end() - 1;
  std::vector<AssignedMaterialSteps>::iterator asIterEnd = assignedSteps.end();

  // loop over the steps
  for (auto& mStep : materialSteps) {
    // start with infinity distance
    double lDist2 = std::numeric_limits<double>::infinity();
    // now assign to a step
    asIter = asIterFlush;
    for (; asIter != asIterEnd; ++asIter) {
      // the currentDist
      double cDist2 = (mStep.position() - asIter->assignedPosition).mag2();
      if (cDist2 < lDist2) {
        // set the new closest distance
        lDist2 = cDist2;
        // remember where you are
        asIterFlush = asIter;
        if (asIterFlush == asIterLast) {
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
          for (auto& sSurface : cLayer->surfaceArray()->surfaces())
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
    if (smp) {
      // get the geo id
      auto   geoID    = surface.geoID();
      size_t volumeID = geoID.value(GeometryID::volume_mask);
      ACTS_VERBOSE("Material surface found with volumeID " << volumeID);
      ACTS_VERBOSE("       - surfaceID is " << geoID.value());
      sMap[geoID] = SurfaceMaterialRecord(surface, *smp->binUtility());
    }
  }
}
