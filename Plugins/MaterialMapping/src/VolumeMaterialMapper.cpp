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

Acts::VolumeMaterialMapper::State
Acts::VolumeMaterialMapper::createState(
    const std::vector<double>& gridAxis1,
    const std::vector<double>& gridAxis2,
    const std::vector<double>& gridAxis3) const
{
  // The Surface material mapping state
  State mState;

  // Add points and sort
  mState.gridPointsPerAxis.push_back(gridAxis1);
  std::sort(mState.gridPointsPerAxis[0].begin(),
            mState.gridPointsPerAxis[0].end());

  // Count the total amount of grid points
  unsigned int numGridPoints = gridAxis1.size();

  // Perform it again if there is a second axis
  if (!gridAxis2.empty()) {
    mState.gridPointsPerAxis.push_back(gridAxis2);
    std::sort(mState.gridPointsPerAxis[1].begin(),
              mState.gridPointsPerAxis[1].end());
    numGridPoints *= gridAxis2.size();

    // Perform it again if there is a third axis
    if (!gridAxis3.empty()) {
      mState.gridPointsPerAxis.push_back(gridAxis3);
      std::sort(mState.gridPointsPerAxis[2].begin(),
                mState.gridPointsPerAxis[2].end());
      numGridPoints *= gridAxis3.size();
    }
  }
  // Set the number of grid points
  mState.accumulatedMaterial.resize(numGridPoints);

  return mState;
}

void
Acts::VolumeMaterialMapper::mapMaterialPoints(
    State&                  mState,
    const RecordedMaterial& mPoints,
    const std::function<unsigned int(const Vector3D&, const State&)>&
        concatenateToGridPoints) const
{
  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    unsigned int index = concatenateToGridPoints(rm.second, mState);
    mState.accumulatedMaterial[index].accumulate(rm.first);
  }
}

std::vector<Acts::Material>
Acts::VolumeMaterialMapper::finalizeMaps(State& mState) const
{
  // Returning material vector
  std::vector<Acts::Material> materialAtGridPoints;
  materialAtGridPoints.reserve(mState.accumulatedMaterial.size());

  // Average material and add it
  for (auto& am : mState.accumulatedMaterial) {
    materialAtGridPoints.push_back(am.totalAverage());
  }
  return materialAtGridPoints;
}