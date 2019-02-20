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
Acts::VolumeMaterialMapper::createState(const vector<double>& gridAxis1,
                                        const vector<double>& gridAxis2,
                                        const vector<double>& gridAxis3) const
{
  // The Surface material mapping state
  State mState;

  // Add points and sort
  mState.gridPointsPerAxis.push_back(gridAxis1);
  std::sort(mState.gridPointsPerAxis[0].begin(),
            mState.gridPointsPerAxis[0].end());

  // Count the total amount of grid points
  unsigned int numGridPoints = gridAxis1.size();

  ACTS_DEBUG("Grid axis found containing " << gridAxis1.size() << " points");

  // Perform it again if there is a second axis
  if (!gridAxis2.empty()) {
    mState.gridPointsPerAxis.push_back(gridAxis2);
    std::sort(mState.gridPointsPerAxis[1].begin(),
              mState.gridPointsPerAxis[1].end());
    numGridPoints *= gridAxis2.size();

    ACTS_DEBUG("Grid axis found containing " << gridAxis2.size() << " points");

    // Perform it again if there is a third axis
    if (!gridAxis3.empty()) {
      mState.gridPointsPerAxis.push_back(gridAxis3);
      std::sort(mState.gridPointsPerAxis[2].begin(),
                mState.gridPointsPerAxis[2].end());
      numGridPoints *= gridAxis3.size();

      ACTS_DEBUG("Grid axis found containing " << gridAxis3.size()
                                               << " points");
    }
  }
  // Set the number of grid points
  mState.accumulatedMaterial.resize(numGridPoints);

  ACTS_DEBUG(mState.accumulatedMaterial.size()
             << " grid points collected ... ");

  return mState;
}

void
Acts::SurfaceMaterialMapper::mapMaterialPoints(
    State&                  mState,
    const RecordedMaterial& mPoints,
    const std::function<unsigned int(const Vector3D&, const State&)>&
        concatenateToGridPoints) const
{
  ACTS_DEBUG("Adding " << mPoints.size() << " material points to the record");

  // Walk over each point
  for (const auto& rm : mPoints) {
    // Search for fitting grid point and accumulate
    unsigned int index = concatenateToGridPoint(rm.second, mState);
    mState.accumulatedMaterial[index].accumulate(rm.first);
  }
}

std::vector<Acts::Material>
Acts::SurfaceMaterialMapper::finalizeMaps(State& mState) const
{
  ACTS_DEBUG("Starting finalization of maps");

  // Returning material vector
  std::vector<Acts::Material> materialAtGridPoint;
  materialAtGridPoint.reserve(mState.accumulatedMaterial.size());

  // Average material and add it
  for (auto& am : mState.accumulatedMaterial) {
    materialAtGridPoint.push_back(am.totalAverage());
  }

  ACTS_VERBOSE("Map finalized");
}