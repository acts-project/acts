// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Digitization/SingleHitSpacePointBuilder.hpp"

template <typename Cluster>
Acts::Vector2D
Acts::SpacePointBuilder<Acts::SingleHitSpacePoint<Cluster>>::localCoords(
    const Cluster& cluster) const {
  // Local position information
  auto par = cluster.parameters();
  Acts::Vector2D local(par[Acts::ParDef::eLOC_0], par[Acts::ParDef::eLOC_1]);
  return local;
}

template <typename Cluster>
Acts::Vector3D
Acts::SpacePointBuilder<Acts::SingleHitSpacePoint<Cluster>>::globalCoords(
    const GeometryContext& gctx, const Cluster& cluster) const {
  // Receive corresponding surface
  auto& clusterSurface = cluster.referenceSurface();

  // Transform local into global position information
  Acts::Vector3D pos, mom;
  clusterSurface.localToGlobal(gctx, localCoords(cluster), mom, pos);

  return pos;
}

template <typename Cluster>
void Acts::SpacePointBuilder<Acts::SingleHitSpacePoint<Cluster>>::
    calculateSpacePoints(const GeometryContext& gctx,
                         const std::vector<const Cluster*>& clusters,
                         std::vector<Acts::SingleHitSpacePoint<Cluster>>&
                             spacePointStorage) const {
  // Set the space point for all stored hits
  for (const auto& c : clusters) {
    Acts::SingleHitSpacePoint<Cluster> spacePoint;
    spacePoint.spacePoint = globalCoords(gctx, *c);
    spacePoint.clusterModule = c;
    spacePointStorage.push_back(std::move(spacePoint));
  }
}