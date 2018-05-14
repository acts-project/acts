// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Tools/SingleHitSpacePointBuilder.hpp"

Acts::Vector2D
Acts::SpacePointBuilder<Acts::SingleHitSpacePoint, void>::localCoords(
    const Acts::PlanarModuleCluster& hit)
{
  // Local position information
  auto           par = hit.parameters();
  Acts::Vector2D local(par[Acts::ParDef::eLOC_0], par[Acts::ParDef::eLOC_1]);
  return local;
}

Acts::Vector3D
Acts::SpacePointBuilder<Acts::SingleHitSpacePoint, void>::globalCoords(
    const Acts::PlanarModuleCluster& hit)
{
  // Receive corresponding surface
  auto& clusterSurface = hit.referenceSurface();

  // Transform local into global position information
  Acts::Vector3D pos, mom;
  clusterSurface.localToGlobal(localCoords(hit), mom, pos);

  return pos;
}

void
Acts::SpacePointBuilder<Acts::SingleHitSpacePoint, void>::addHits(
    std::vector<Acts::SingleHitSpacePoint>&              spacePointStorage,
    const std::vector<Acts::PlanarModuleCluster const*>& hits)
{
  // Walk over every hit and add them
  for (auto& cluster : hits) {
    // Declare helper variable
    Acts::SingleHitSpacePoint tmpSpacePoint;
    tmpSpacePoint.hitModule.push_back(cluster);
    spacePointStorage.push_back(tmpSpacePoint);
  }
}

void
Acts::SpacePointBuilder<Acts::SingleHitSpacePoint, void>::calculateSpacePoints(
    std::vector<Acts::SingleHitSpacePoint>& spacePointStorage)
{
  // Set the space point for all stored hits
  for (auto& sPoint : spacePointStorage) {
    if (sPoint.spacePoint != Acts::Vector3D::Zero(3)) continue;
    if (sPoint.hitModule.size() != 1) continue;
    sPoint.spacePoint = globalCoords(*(sPoint.hitModule[0]));
  }
}
