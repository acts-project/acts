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
Acts::OneHitSpacePointBuilder::addHits(
    const std::vector<Acts::PlanarModuleCluster>& vec)
{
  // Declare helper variable
  std::shared_ptr<Acts::SpacePoint> tmpSpacePoint;

  // Add hits
  for (auto& cluster : vec) {
    tmpSpacePoint->hitModule = &cluster;
    m_allSpacePoints.push_back(tmpSpacePoint);
  }
}

void
Acts::OneHitSpacePointBuilder::addSpacePoint(
    const std::shared_ptr<Acts::SpacePoint>& sPoint)
{
  // Add a hit if the module is set
  if (sPoint->hitModule) m_allSpacePoints.push_back(sPoint);
}

void
Acts::OneHitSpacePointBuilder::calculateSpacePoints()
{
  // Set the space point for all stored hits
  for (auto& sPoint : m_allSpacePoints) {
    if (sPoint->spacePoint != Acts::Vector3D::Zero(3)) continue;
    sPoint->spacePoint = globalCoords(*(sPoint->hitModule));
  }
}

const std::vector<std::shared_ptr<Acts::SpacePoint>>&
Acts::OneHitSpacePointBuilder::spacePoints()
{
  return m_allSpacePoints;
}
