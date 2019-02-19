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

#include "Acts/Plugins/MaterialMapping/VolumeMaterialMapper.hpp"

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
  std::sort(mState.edgesPerAxis[0].begin(), mState.edgesPerAxis[0].end());
  if(!edgeAxis2.empty())
  {
	mState.edgesPerAxis.push_back(edgeAxis2);
	std::sort(mState.edgesPerAxis[1].begin(), mState.edgesPerAxis[1].end());
	numEdges += edgeAxis2.size();
	}
  if(!edgeAxis3.empty())
  {
	mState.edgesPerAxis.push_back(edgeAxis3);
	std::sort(mState.edgesPerAxis[2].begin(), mState.edgesPerAxis[2].end());
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

std::vector<Acts::Material>
Acts::SurfaceMaterialMapper::finalizeMaps(State& mState) const
{
	std::vector<Acts::Material> materialAtEdge;
	materialAtEdge.reserve(mState.accumulatedMaterial.size());

  for (auto& am : mState.accumulatedMaterial) {
	  materialAtEdge.push_back(am.totalAverage().first.material());
  }
}

// TODO: note down that this function treats everything as if it would be ordered
void
Acts::SurfaceMaterialMapper::mapMaterialTrack(
    State&                       mState,
    const RecordedMaterialTrack& mTrack, const std::function<unsigned int(const Vector3D&, const State&)>& concatenateToEdge) const
{
  for(const auto& rmp : mTrack)
  {
	unsigned int index = concatenateToEdge(rmp.second, mState);
	mState.accumulatedMaterial[index].accumulate(rmp.second, rmp.first);
  }
  
  for(auto& am : mState.accumulatedMaterial)
  {
	  am.eventAverage();
  }
}
