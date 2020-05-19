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

  for (auto& sf : surfaces) {
    auto sfExtent = sf->polyhedronRepresentation(gctx, 1).extent();
  }

  return protoLayers;
}