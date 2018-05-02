// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Tools/OneHitSpacePointBuilder.hpp"

Acts::OneHitSpacePointBuilder::OneHitSpacePointBuilder()
{
}

Acts::Vector2D
Acts::OneHitSpacePointBuilder::localCoords(
    const Acts::PlanarModuleCluster& hit) const
{
  // Local position information
  auto           par = hit.parameters();
  Acts::Vector2D local(par[Acts::ParDef::eLOC_0], par[Acts::ParDef::eLOC_1]);
  return local;
}

Acts::Vector3D
Acts::OneHitSpacePointBuilder::globalCoords(
    const Acts::PlanarModuleCluster& hit) const
{
  // Receive corresponding surface
  auto& clusterSurface = hit.referenceSurface();

  // Transform local into global position information
  Acts::Vector3D pos, mom;
  clusterSurface.localToGlobal(localCoords(hit), mom, pos);

  return pos;
}

void
Acts::OneHitSpacePointBuilder::addHits(std::vector<std::vector<Acts::PlanarModuleCluster const*>>& hits)
{
  // Check that only a single surface is used
  assert(hits.size() == 1);
  
  // Walk over every hit and add them
  for (auto& cluster : hits[0]) {
	// Declare helper variable
    Acts::SpacePoint tmpSpacePoint;
	tmpSpacePoint.hitModule.resize(1);
	tmpSpacePoint.hitModule.push_back(cluster);
    m_allSpacePoints.push_back(tmpSpacePoint);
  }
}

void
Acts::OneHitSpacePointBuilder::addSpacePoint(
    Acts::SpacePoint& sPoint)
{
  // Add a hit if the module is set
  if (sPoint.hitModule.size() == 1 && sPoint.hitModule[0]) m_allSpacePoints.push_back(sPoint);
}

void
Acts::OneHitSpacePointBuilder::calculateSpacePoints()
{
  // Set the space point for all stored hits
  for (auto& sPoint : m_allSpacePoints) {
    if (sPoint.spacePoint != Acts::Vector3D::Zero(3)) continue;
    sPoint.spacePoint = globalCoords(*(sPoint.hitModule[0]));
  }
}

const std::vector<Acts::SpacePoint>&
Acts::OneHitSpacePointBuilder::spacePoints()
{
  return m_allSpacePoints;
}
