// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/VolumeMaterialMapper.hpp"
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/MaterialMapUtils.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Material/InterpolatedMaterialMap.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Utilities/BinAdjustmentVolume.hpp"

namespace {
using EAxis = Acts::detail::EquidistantAxis;
using Grid2D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis>;
using Grid3D =
    Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
using MaterialGrid2D = Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis>;
using MaterialGrid3D =
    Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis, EAxis>;

}  // namespace

Acts::VolumeMaterialMapper::VolumeMaterialMapper(
    const Config& cfg, StraightLinePropagator propagator,
     std::unique_ptr<const Logger> slogger)
    : m_cfg(cfg),
      m_propagator(std::move(propagator)),
      m_logger(std::move(slogger)) {}

Acts::VolumeMaterialMapper::State Acts::VolumeMaterialMapper::createState(
    const GeometryContext& gctx, const MagneticFieldContext& mctx,
    const TrackingGeometry& tGeometry) const {
  // Parse the geometry and find all surfaces with material proxies
  auto world = tGeometry.highestTrackingVolume();

  // The Surface material mapping state
  State mState(gctx, mctx);
  resolveMaterialVolume(mState, *world);
  collectMaterialSurface(mState, *world);

  return mState;
}

void Acts::VolumeMaterialMapper::resolveMaterialVolume(
    State& mState, const TrackingVolume& tVolume) const {
  ACTS_VERBOSE("Checking volume '" << tVolume.volumeName()
                                   << "' for material surfaces.")

  ACTS_VERBOSE("- Insert Volume ...");
  checkAndInsert(mState, tVolume);

  // Step down into the sub volume
  if (tVolume.confinedVolumes()) {
    ACTS_VERBOSE("- Check children volume ...");
    for (auto& sVolume : tVolume.confinedVolumes()->arrayObjects()) {
      // Recursive call
      resolveMaterialVolume(mState, *sVolume);
    }
  }
  if (!tVolume.denseVolumes().empty()) {
    for (auto& sVolume : tVolume.denseVolumes()) {
      // Recursive call
      resolveMaterialVolume(mState, *sVolume);
    }
  }
}

void Acts::VolumeMaterialMapper::checkAndInsert(State& mState,
                                                 const TrackingVolume& volume) const
                                                 {
  auto volumeMaterial = volume.volumeMaterial();
  // Check if the volume has a proxy
  if (volumeMaterial != nullptr) {

    auto geoID = volume.geoID();
    size_t volumeID = geoID.volume();
    ACTS_DEBUG("Material volume found with volumeID " << volumeID);
    ACTS_DEBUG("       - ID is " << geoID);

    RecordedMaterialPoint mat;
    mState.recordedMaterial[geoID] = mat;

    // We need a dynamic_cast to either a volume material proxy or
    // proper surface material
    auto psm = dynamic_cast<const ProtoVolumeMaterial*>(volumeMaterial);
    // Get the bin utility: try proxy material first
    const BinUtility* bu = (psm != nullptr) ? (&psm->binUtility()) : nullptr;
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - (proto) binning is " << *bu);
      // Now update
      BinUtility buAdjusted = adjustBinUtility(*bu, volume);
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - adjusted binning is " << buAdjusted);
      mState.materialBin[geoID] = buAdjusted;
      return;
    }
    //TODO!!!
    // Second attempt: binned material
    auto bmp = dynamic_cast<const InterpolatedMaterialMap<MaterialMapper<MaterialGrid3D>>*>(volumeMaterial);
    bu = (bmp != nullptr) ? (&bmp->binUtility()) : nullptr;
    // Creaete a binned type of material
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_DEBUG("       - binning is " << *bu);
      mState.materialBin[geoID] = *bu;
      return;
    } else {
      // Create a homogeneous type of material
      ACTS_DEBUG("       - this is homogeneous material.");
      return;
    }
  }
}

void Acts::VolumeMaterialMapper::collectMaterialSurface(
    State& mState, const TrackingVolume& tVolume) const {
  ACTS_VERBOSE("Checking volume '" << tVolume.volumeName()
                                   << "' for material surfaces.")

  ACTS_VERBOSE("- boundary surfaces ...");
  // Check the boundary surfaces
  for (auto& bSurface : tVolume.boundarySurfaces()) {
    if(bSurface->surfaceRepresentation().surfaceMaterial() != nullptr){
      mState.surfaceMaterial[bSurface->surfaceRepresentation().geoID()] = bSurface->surfaceRepresentation().surfaceMaterialSharedPtr();
    }
  }

  ACTS_VERBOSE("- confined layers ...");
  // Check the confined layers
  if (tVolume.confinedLayers() != nullptr) {
    for (auto& cLayer : tVolume.confinedLayers()->arrayObjects()) {
      // Take only layers that are not navigation layers
      if (cLayer->layerType() != navigation) {
        // Check the representing surface
        if(cLayer->surfaceRepresentation().surfaceMaterial() != nullptr){
          mState.surfaceMaterial[cLayer->surfaceRepresentation().geoID()] = cLayer->surfaceRepresentation().surfaceMaterialSharedPtr();
        }
        // Get the approach surfaces if present
        if (cLayer->approachDescriptor() != nullptr) {
          for (auto& aSurface :
               cLayer->approachDescriptor()->containedSurfaces()) {
            if (aSurface != nullptr) {
              if(aSurface->surfaceMaterial() != nullptr){
                mState.surfaceMaterial[aSurface->geoID()] = aSurface->surfaceMaterialSharedPtr();
              }
            }
          }
        }
        // Get the sensitive surface is present
        if (cLayer->surfaceArray() != nullptr) {
          // Sensitive surface loop
          for (auto& sSurface : cLayer->surfaceArray()->surfaces()) {
            if (sSurface != nullptr) {
              if(sSurface->surfaceMaterial() != nullptr){
                mState.surfaceMaterial[sSurface->geoID()] = sSurface->surfaceMaterialSharedPtr();
              }
            }
          }
        }
      }
    }
  }
  // Step down into the sub volume
  if (tVolume.confinedVolumes()) {
    for (auto& sVolume : tVolume.confinedVolumes()->arrayObjects()) {
      // Recursive call
      collectMaterialSurface(mState, *sVolume);
    }
  }
}

void Acts::VolumeMaterialMapper::finalizeMaps(State& mState) const {
  // iterate over the map to call the total average
  for (auto& recMaterial : mState.recordedMaterial) {

    if(mState.materialBin[recMaterial.first].dimensions() == 0){
      Acts::AccumulatedVolumeMaterial homogeneousAccumulation;
      for (const auto& rm : recMaterial.second) {
        homogeneousAccumulation.accumulate(rm.first);
      }
      Acts::Material mat = homogeneousAccumulation.average();
      mState.volumeMaterial[recMaterial.first] = std::make_unique<HomogeneousVolumeMaterial>(std::move(mat));
      return;
    }
    else if (mState.materialBin[recMaterial.first].dimensions() == 3){
      ACTS_DEBUG("Create the material for volume  " << recMaterial.first);
      std::function<Acts::Vector3D(Acts::Vector3D)> transfoGlobalToLocal;
      std::function<Grid3D::index_t(const Acts::Vector3D&, const Grid3D&)> mapMaterial;
      Grid3D Grid = createGrid(mState.materialBin[recMaterial.first], transfoGlobalToLocal, mapMaterial);
      MaterialGrid3D matGrid = mapMaterialPoints(Grid, recMaterial.second, mapMaterial);
      MaterialMapper<MaterialGrid3D> matMap(transfoGlobalToLocal, matGrid);
      mState.volumeMaterial[recMaterial.first] = std::make_unique<InterpolatedMaterialMap<MaterialMapper<MaterialGrid3D>>>(std::move(matMap),std::move(mState.materialBin[recMaterial.first]));
      return;
    }
    else{
      throw std::invalid_argument("Incorrect bin dimension, only 0 and 3 are accepted");
    }
  }
}

void Acts::VolumeMaterialMapper::mapMaterialTrack(
    State& mState, RecordedMaterialTrack& mTrack) const {
  // Neutral curvilinear parameters
  NeutralCurvilinearParameters start(std::nullopt, mTrack.first.first,
                                     mTrack.first.second, 0.);

  // Prepare Action list and abort list
  using DebugOutput = detail::DebugOutputActor;
  using MaterialVolumeCollector = VolumeCollector<MaterialVolume>;
  using ActionList = ActionList<MaterialVolumeCollector, DebugOutput>;
  using AbortList = AbortList<detail::EndOfWorldReached>;

  PropagatorOptions<ActionList, AbortList> options(mState.geoContext,
                                                   mState.magFieldContext);
  options.debug = m_cfg.mapperDebugOutput;

  // Now collect the material layers by using the straight line propagator
  const auto& result = m_propagator.propagate(start, options).value();
  auto mcResult = result.get<MaterialVolumeCollector::result_type>();
  // Massive screen output
  if (m_cfg.mapperDebugOutput) {
    auto debugOutput = result.get<DebugOutput::result_type>();
    ACTS_VERBOSE("Debug propagation output.");
    ACTS_VERBOSE(debugOutput.debugString);
  }

  auto mappingVolumes = mcResult.collected;

  // Retrieve the recorded material from the recorded material track
  auto& rMaterial = mTrack.second.materialInteractions;
  ACTS_VERBOSE("Retrieved " << rMaterial.size()
                            << " recorded material steps to map.")

  // These should be mapped onto the mapping surfaces found
  ACTS_VERBOSE("Found     " << mappingVolumes.size()
                            << " mapping volumes for this track.");
  ACTS_VERBOSE("Mapping volumes are :")
  for (auto& mVolumes : mappingVolumes) {
    ACTS_VERBOSE(" - Volume : " << mVolumes.volume->geoID()
                                 << " at position = (" <<
                                 mVolumes.position.x()
                                 << ", " << mVolumes.position.y() << ", "
                                 << mVolumes.position.z() << ")");

    mappingVolumes.push_back(mVolumes);
  }
  auto rmIter = rMaterial.begin();
  auto volIter = mappingVolumes.begin();
  bool encounterVolume = false;
  double mappingStep = 1.;

  int volumeStep = 1;
  Acts::Vector3D extraPosition = {0,0,0};
  Acts::Vector3D extraDirection = {0,0,0};


  while (rmIter != rMaterial.end() && volIter != mappingVolumes.end()) {
    if(volIter != mappingVolumes.end() && encounterVolume==true && !volIter->volume->inside(rmIter->position)){
      encounterVolume=false;
      ++volIter;
    }

    if(volIter != mappingVolumes.end() && volIter->volume->inside(rmIter->position)){
      volumeStep = floor(rmIter->materialProperties.thickness()/mappingStep);
      auto properties = rmIter->materialProperties;
      float remainder = rmIter->materialProperties.thickness() - mappingStep*volumeStep;
      properties.scaleThickness(mappingStep/properties.thickness());
      mState.recordedMaterial[volIter->volume->geoID()].push_back(std::pair(properties,rmIter->position));
      for(int step=1;step<=volumeStep;step++){
        extraDirection = rmIter->direction;
        extraDirection = extraDirection*(mappingStep/extraDirection.norm());
        extraPosition = rmIter->position + step*extraDirection;
        if(step==volumeStep){
          properties.scaleThickness(remainder/properties.thickness());
        }
        mState.recordedMaterial[volIter->volume->geoID()].push_back(std::pair(properties,extraPosition));
      }
      encounterVolume=true;
    }
    ++rmIter;
  }
}
