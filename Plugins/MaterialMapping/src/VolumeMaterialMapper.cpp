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
    std::unique_ptr<const Logger> slogger)
  : m_logger(std::move(slogger))
{
}

Acts::VolumeMaterialMapper::State
Acts::VolumeMaterialMapper::createState(
    const vector<double>& gridAxis1, const vector<double>& gridAxis2, const vector<double>& gridAxis3) const
{
  // The Surface material mapping state
  State mState;
  
  unsigned int numGridPoints = gridAxis1.size();
  mState.gridPointsPerAxis.push_back(gridAxis1);
  std::sort(mState.gridPointsPerAxis[0].begin(), mState.gridPointsPerAxis[0].end());
  if(!gridAxis2.empty())
  {
	mState.gridPointsPerAxis.push_back(gridAxis2);
	std::sort(mState.gridPointsPerAxis[1].begin(), mState.gridPointsPerAxis[1].end());
	numGridPoints *= gridAxis2.size();
	}
  if(!gridAxis3.empty())
  {
	mState.gridPointsPerAxis.push_back(gridAxis3);
	std::sort(mState.gridPointsPerAxis[2].begin(), mState.gridPointsPerAxis[2].end());
	numGridPoints *= gridAxis3.size();
	}

	mState.accumulatedMaterial.resize(numGridPoints);

  ACTS_DEBUG(mState.accumulatedMaterial.size()
             << " anchor points collected ... ");
  ACTS_VERBOSE(numGridPoints << " anchor points found");
 
  return mState;
}

std::vector<Acts::Material>
Acts::SurfaceMaterialMapper::finalizeMaps(State& mState) const
{
	std::vector<Acts::Material> materialAtGridPoint;
	materialAtGridPoint.reserve(mState.accumulatedMaterial.size());

  for (auto& am : mState.accumulatedMaterial) {
	  materialAtGridPoint.push_back(am.totalAverage().first.material());
  }
}

void
Acts::SurfaceMaterialMapper::mapMaterialPoints(
    State&                       mState,
    const RecordedMaterialTrack& mPoints, const std::function<unsigned int(const Vector3D&, const State&)>& concatenateToGridPoints) const
{
  for(const auto& rmp : mPoints)
  {
	unsigned int index = concatenateToGridPoint(rmp.second, mState);
	mState.accumulatedMaterial[index].accumulate(rmp.second, rmp.first);
  }
}
