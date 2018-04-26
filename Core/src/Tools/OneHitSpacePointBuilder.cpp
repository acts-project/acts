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
  Acts::OneHitSpacePointBuilder::Hit tmpHits;

  // Add hits
  for (auto& cluster : vec) {
    tmpHits.hitModule = &cluster;
    m_allHits.push_back(tmpHits);
  }
}

void
Acts::OneHitSpacePointBuilder::addHit(
    const Acts::OneHitSpacePointBuilder::Hit& hit)
{
  // Add a hit if the module is set
  if (hit.hitModule) m_allHits.push_back(hit);
}

void
Acts::OneHitSpacePointBuilder::calculateSpacePoints()
{
  // Set the space point for all stored hits
  for (auto& hit : m_allHits) {
    if (hit.spacePoint != Acts::Vector3D::Zero(3)) continue;
    hit.spacePoint = globalCoords(*(hit.hitModule));
  }
}

const std::vector<Acts::OneHitSpacePointBuilder::Hit>&
Acts::OneHitSpacePointBuilder::hits()
{
  return m_allHits;
}
