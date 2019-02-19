// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// VolumeMaterialMapper.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Plugins/MaterialMapping/SurfaceMaterialMapper.hpp"
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Material/SurfaceMaterialProxy.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"

Acts::VolumeMaterialMapper::VolumeMaterialMapper(
    bool mapperDebugOutput                 dbg,
    std::unique_ptr<const Logger> slogger)
  : m_mapperDebugOutput(dpg)
  , m_logger(std::move(slogger))
{
}

Acts::VolumeMaterialMapper::State
Acts::VolumeMaterialMapper::createState(
    const vector<double>& edgeAxis1, const vector<double>& edgeAxis2, const vector<double>& edgeAxis3) const
{
  // The Surface material mapping state
  State mState;
  
  unsigned int numEdges = edgeAxis1.size();
  mState.edgesPerAxis.push_back(edgeAxis1);
  if(!edgeAxis2.empty())
  {
	mState.edgesPerAxis.push_back(edgeAxis2);
	numEdges += edgeAxis2.size();
	}
  if(!edgeAxis3.empty())
  {
	mState.edgesPerAxis.push_back(edgeAxis3);
	numEdges += edgeAxis3.size();
	}

	mState.accumulatedMaterial.resize(numEdges);

  ACTS_DEBUG(mState.accumulatedMaterial.size()
             << " anchor points collected ... ");
  for (auto& axis : mState.edgesPerAxis) {
    ACTS_VERBOSE(" -> Anchor points per axis: " << axis.size());
  }
  return mState;
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
      nullptr, mTrack.position(), mTrack.direction());

  // Prepare Action list and abort list
  using DebugOutput              = detail::DebugOutputActor;
  using MaterialSurfaceCollector = SurfaceCollector<MaterialSurface>;
  using ActionList = ActionList<MaterialSurfaceCollector, DebugOutput>;
  using AbortList  = AbortList<detail::EndOfWorldReached>;

  PropagatorOptions<ActionList, AbortList> options;
  options.debug = m_mapperDebugOutput;

  // Now collect the material layers by using the straight line propagator
  const auto& result   = m_propagator.propagate(start, options);
  auto        mcResult = result.get<MaterialSurfaceCollector::result_type>();
  auto        mappingSurfaces = mcResult.collected;

  // Retrieve the recorded material
  const auto& rMaterial = mTrack.recordedMaterialProperties();

  ACTS_VERBOSE("Retrieved " << rMaterial.size()
                            << " recorded material properties to map.")

  ACTS_VERBOSE("Found     " << mappingSurfaces.size()
                            << " mapping surfaces for this track.");

  // Prepare the assignment store
  std::vector<AssignedMaterialProperties> assignedMaterial;
  assignedMaterial.reserve(mappingSurfaces.size());
  for (auto& mSurface : mappingSurfaces) {
    Intersection msIntersection = mSurface.surface->intersectionEstimate(
        mTrack.position(), mTrack.direction(), forward, true);
    if (msIntersection) {
      double pathCorrection = mSurface.surface->pathCorrection(
          msIntersection.position, mTrack.direction());
      AssignedMaterialProperties amp(
          mSurface.surface->geoID(), msIntersection.position, pathCorrection);
      assignedMaterial.push_back(std::move(amp));
    }
  }

  ACTS_VERBOSE("Prepared  " << assignedMaterial.size()
                            << " assignment stores for this event.");

  if (!assignedMaterial.empty()) {
    // Match the recorded material to the assigment stores
    auto aStore = assignedMaterial.begin();
    // This assumes ordered recorded material
    for (auto rmp : rMaterial) {
      // if it's not the last & the next one is closer : switch
      if (aStore != assignedMaterial.end() - 1) {
        if ((aStore->assignedPosition - rmp.second).norm()
            > ((aStore + 1)->assignedPosition - rmp.second).norm()) {
          ++aStore;
        }
      }
      // Now assign
      aStore->assignedProperties.push_back(rmp);
    }

    // Now move the assigned properties into accumulation map
    for (auto aprop : assignedMaterial) {
      /// get the according map
      auto aSurfaceMaterial = mState.accumulatedMaterial.find(aprop.geoID);
      // you have assigned material
      if (!aprop.assignedProperties.empty()) {
        aSurfaceMaterial->second.accumulate(aprop.assignedPosition,
                                            aprop.assignedProperties,
                                            1. / aprop.pathCorrection);
      } else {
        // assign a single vacuum step to regulate the (correct) averaging
        aSurfaceMaterial->second.accumulate(aprop.assignedPosition,
                                            MaterialProperties{1.});
      }
      // now average over the event
      aSurfaceMaterial->second.eventAverage();
    }
  }
}
