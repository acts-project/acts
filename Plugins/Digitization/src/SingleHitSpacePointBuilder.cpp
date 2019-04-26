// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Digitization/SingleHitSpacePointBuilder.hpp"

Acts::Vector2D Acts::SpacePointBuilder<Acts::SingleHitSpacePoint>::localCoords(
    const PlanarModuleCluster& cluster) const {
  // Local position information
  auto par = cluster.parameters();
  Acts::Vector2D local(par[Acts::ParDef::eLOC_0], par[Acts::ParDef::eLOC_1]);
  return local;
}

Acts::Vector3D Acts::SpacePointBuilder<Acts::SingleHitSpacePoint>::globalCoords(
    const GeometryContext& gctx, const PlanarModuleCluster& cluster) const {
  // Receive corresponding surface
  auto& clusterSurface = cluster.referenceSurface();

  // Transform local into global position information
  Acts::Vector3D pos, mom;
  clusterSurface.localToGlobal(gctx, localCoords(cluster), mom, pos);

  return pos;
}

void Acts::SpacePointBuilder<Acts::SingleHitSpacePoint>::calculateSpacePoints(
    const GeometryContext& gctx,
    const std::vector<const Acts::PlanarModuleCluster*>& clusters,
    std::vector<Acts::SingleHitSpacePoint>& spacePointStorage) const {
  // Set the space point for all stored hits
  for (const auto& c : clusters) {
    Acts::SingleHitSpacePoint spacePoint;
    spacePoint.spacePoint = globalCoords(gctx, *c);
    spacePoint.clusterModule = c;
    spacePointStorage.push_back(std::move(spacePoint));
  }
}
