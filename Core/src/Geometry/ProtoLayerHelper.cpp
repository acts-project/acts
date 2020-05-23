// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/ProtoLayerHelper.hpp"
#include "Acts/Geometry/Extent.hpp"

std::vector<Acts::ProtoLayer> Acts::ProtoLayerHelper::protoLayers(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    BinningValue bValue, double joinTolerance) const {
  std::vector<Acts::ProtoLayer> protoLayers;

  using SurfaceCluster = std::pair<Extent, std::vector<const Surface*>>;
  std::vector<SurfaceCluster> clusteredSurfaces;

  /// Helper function to find/create the cluster of surfaces where
  /// the Extent belongs to. In case none is found, a new one is inserted
  ///
  /// @param extent The test extent for finding the cluster
  ///
  /// @return the referece of the SurfaceCluster for insertion
  auto findCluster = [&](const Extent& extent) -> SurfaceCluster& {
    for (auto& cluster : clusteredSurfaces) {
      if (cluster.first.intersects(extent, bValue, joinTolerance)) {
        return cluster;
      }
    }
    // No cluster found, let's create a new one
    clusteredSurfaces.push_back(SurfaceCluster(extent, {}));
    return clusteredSurfaces.back();
  };

  // Loop over surfaces and sort into clusters
  for (auto& sf : surfaces) {
    auto sfExtent = sf->polyhedronRepresentation(gctx, 1).extent();
    auto sfCluster = findCluster(sfExtent);
    sfCluster.first.extend(sfExtent);
    sfCluster.second.push_back(sf);
  }

  // Loop over clusters and create ProtoLayer
  protoLayers.reserve(clusteredSurfaces.size());
  for (auto& clusters : clusteredSurfaces) {
    protoLayers.push_back(ProtoLayer(gctx, clusters.second));
  }

  return protoLayers;
}