// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialMapper.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Plugins/MaterialMapping/SurfaceMaterialMapper.hpp"
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Utilities/BinAdjustment.hpp"
#include "Acts/Utilities/BinUtility.hpp"

Acts::SurfaceMaterialMapper::SurfaceMaterialMapper(
    const Config&                 cfg,
    StraightLinePropagator        propagator,
    std::unique_ptr<const Logger> slogger)
  : m_cfg(cfg)
  , m_propagator(std::move(propagator))
  , m_logger(std::move(slogger))
{
}

Acts::SurfaceMaterialMapper::State
Acts::SurfaceMaterialMapper::createState(
    const GeometryContext&      gctx,
    const MagneticFieldContext& mctx,
    const TrackingGeometry&     tGeometry) const
{
  // Parse the geometry and find all surfaces with material proxies
  auto world = tGeometry.highestTrackingVolume();

  // The Surface material mapping state
  State mState(gctx, mctx);
  resolveMaterialSurfaces(mState, *world);

  ACTS_DEBUG(mState.accumulatedMaterial.size()
             << " Surfaces with PROXIES collected ... ");
  for (auto& smg : mState.accumulatedMaterial) {
    ACTS_VERBOSE(" -> Surface in with id " << smg.first.toString());
  }
  return mState;
}

void
Acts::SurfaceMaterialMapper::resolveMaterialSurfaces(
    State&                mState,
    const TrackingVolume& tVolume) const
{
  ACTS_VERBOSE("Checking volume '" << tVolume.volumeName()
                                   << "' for material surfaces.")
  // check the boundary surfaces
  for (auto& bSurface : tVolume.boundarySurfaces()) {
    checkAndInsert(mState, bSurface->surfaceRepresentation());
  }
  // check the confined layers
  if (tVolume.confinedLayers() != nullptr) {
    for (auto& cLayer : tVolume.confinedLayers()->arrayObjects()) {
      // take only layers that are not navigation layers
      if (cLayer->layerType() != navigation) {
        // check the representing surface
        checkAndInsert(mState, cLayer->surfaceRepresentation());
        // get the approach surfaces if present
        if (cLayer->approachDescriptor() != nullptr) {
          for (auto& aSurface :
               cLayer->approachDescriptor()->containedSurfaces()) {
            if (aSurface != nullptr) {
              checkAndInsert(mState, *aSurface);
            }
          }
        }
        // get the sensitive surface is present
        if (cLayer->surfaceArray() != nullptr) {
          // sensitive surface loop
          for (auto& sSurface : cLayer->surfaceArray()->surfaces()) {
            if (sSurface != nullptr) {
              checkAndInsert(mState, *sSurface);
            }
          }
        }
      }
    }
  }
  // step down into the sub volume
  if (tVolume.confinedVolumes()) {
    for (auto& sVolume : tVolume.confinedVolumes()->arrayObjects()) {
      // recursive call
      resolveMaterialSurfaces(mState, *sVolume);
    }
  }
}

void
Acts::SurfaceMaterialMapper::checkAndInsert(State&         mState,
                                            const Surface& surface) const
{

  auto surfaceMaterial = surface.surfaceMaterial();
  // check if the surface has a proxy
  if (surfaceMaterial != nullptr) {

    auto   geoID    = surface.geoID();
    size_t volumeID = geoID.value(GeometryID::volume_mask);
    ACTS_INFO("Material surface found with volumeID " << volumeID);
    ACTS_INFO("       - surfaceID is " << geoID.toString());

    // We need a dynamic_cast to either a surface material proxy or
    // proper surface material
    auto psm = dynamic_cast<const ProtoSurfaceMaterial*>(surfaceMaterial);

    // Get the bin utility: try proxy material first
    const BinUtility* bu = (psm != nullptr) ? (&psm->binUtility()) : nullptr;
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_INFO("       - (proto) binning is " << *bu);
      // Now update
      BinUtility buAdjusted = adjustBinUtility(*bu, surface);
      // Screen output for Binned Surface material
      ACTS_INFO("       - adjusted binning is " << buAdjusted);
      mState.accumulatedMaterial[geoID]
          = AccumulatedSurfaceMaterial(buAdjusted);
      return;
    }

    // Second attempt: binned material
    auto bmp = dynamic_cast<const BinnedSurfaceMaterial*>(surfaceMaterial);
    bu       = (bmp != nullptr) ? (&bmp->binUtility()) : nullptr;
    // Creaete a binned type of material
    if (bu != nullptr) {
      // Screen output for Binned Surface material
      ACTS_INFO("       - binning is " << *bu);
      mState.accumulatedMaterial[geoID] = AccumulatedSurfaceMaterial(*bu);
    } else {
      // Create a homogeneous type of material
      ACTS_INFO("       - this is homogeneous material.");
      mState.accumulatedMaterial[geoID] = AccumulatedSurfaceMaterial();
    }
  }
}

void
Acts::SurfaceMaterialMapper::finalizeMaps(State& mState) const
{
  // iterate over the map to call the total average
  for (auto& accMaterial : mState.accumulatedMaterial) {
    mState.surfaceMaterial[accMaterial.first]
        = accMaterial.second.totalAverage();
  }
}

void
Acts::SurfaceMaterialMapper::mapMaterialTrack(
    State&                       mState,
    const RecordedMaterialTrack& mTrack) const
{
  // Neutral curvilinear parameters
  NeutralCurvilinearParameters start(
      nullptr, mTrack.first.first, mTrack.first.second);

  // Prepare Action list and abort list
  using DebugOutput              = detail::DebugOutputActor;
  using MaterialSurfaceCollector = SurfaceCollector<MaterialSurface>;
  using ActionList = ActionList<MaterialSurfaceCollector, DebugOutput>;
  using AbortList  = AbortList<detail::EndOfWorldReached>;

  PropagatorOptions<ActionList, AbortList> options(mState.geoContext,
                                                   mState.magFieldContext);
  options.debug = m_cfg.mapperDebugOutput;

  // Now collect the material layers by using the straight line propagator
  const auto& result   = m_propagator.propagate(start, options).value();
  auto        mcResult = result.get<MaterialSurfaceCollector::result_type>();
  auto        mappingSurfaces = mcResult.collected;

  // Retrieve the recorded material from the recorded material track
  const auto& rMaterial = mTrack.second.materialInteractions;
  ACTS_VERBOSE("Retrieved " << rMaterial.size()
                            << " recorded material properties to map.")

  // These should be mapped onto the mapping surfaces found
  ACTS_VERBOSE("Found     " << mappingSurfaces.size()
                            << " mapping surfaces for this track.");

  // Run the mapping process, i.e. take the recorded material and map it
  // onto the mapping surfaces
  //
  // The material steps and surfaces are assumed to be ordered along the
  // mapping ray:
  auto rmIter = rMaterial.begin();
  auto sfIter = mappingSurfaces.begin();

  // Use those to minimize the lookup
  GeometryID lastID    = GeometryID();
  GeometryID currentID = GeometryID();
  Vector3D   currentPos(0., 0., 0);
  double     currentPathCorrection = 0.;
  auto       currentAccMaterial    = mState.accumulatedMaterial.end();

  // To remember the bins of this event
  using MapBin = std::pair<AccumulatedSurfaceMaterial*, std::array<size_t, 3>>;
  std::multimap<AccumulatedSurfaceMaterial*, std::array<size_t, 3>>
      touchedMapBins;

  // Assign the recorded ones, break if you hit an end
  while (rmIter != rMaterial.end() && sfIter != mappingSurfaces.end()) {
    // First check if the distance to the next surface is already closer
    // don't do the check for the last one, stay on the last possible surface
    if (sfIter != mappingSurfaces.end() - 1
        && (rmIter->position - sfIter->position).norm()
            > (rmIter->position - (sfIter + 1)->position).norm()) {
      // switch to next assignment surface
      // @TODO: empty hits, i.e. surface is hit but,
      // has no recorded material assigned
      ++sfIter;
    }
    // get the current Surface ID
    currentID = sfIter->surface->geoID();
    // We have work to do: the assignemnt surface has changed
    if (currentID != lastID) {
      // Let's (re-)assess the information
      lastID                = currentID;
      currentPos            = (sfIter)->position;
      currentPathCorrection = sfIter->surface->pathCorrection(
          mState.geoContext, currentPos, sfIter->direction);
      currentAccMaterial = mState.accumulatedMaterial.find(currentID);
    }
    // Now assign the material for the accumulation process
    auto tBin = currentAccMaterial->second.accumulate(
        currentPos, rmIter->materialProperties, currentPathCorrection);
    touchedMapBins.insert(MapBin(&(currentAccMaterial->second), tBin));
    // Switch to next material
    ++rmIter;
  }

  // After mapping this track, average the touched bins
  for (auto tmapBin : touchedMapBins) {
    tmapBin.first->trackAverage({tmapBin.second});
  }
}
