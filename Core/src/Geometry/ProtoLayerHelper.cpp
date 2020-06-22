// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <map>

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/ProtoLayerHelper.hpp"

std::vector<Acts::ProtoLayer> Acts::ProtoLayerHelper::protoLayers(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    const SortingConfig& sorting,
    const std::vector<BinningConfig>& binning) const {
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
      if (cluster.first.intersects(extent, sorting.first, sorting.second)) {
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
    auto& sfCluster = findCluster(sfExtent);
    sfCluster.first.extend(sfExtent);
    sfCluster.second.push_back(sf);
  }

  // Loop over clusters and create ProtoLayer
  protoLayers.reserve(clusteredSurfaces.size());
  for (auto& clusters : clusteredSurfaces) {
    ACTS_VERBOSE("Creating ProtoLayer with " << clusters.second.size()
                                             << " surfaces.");
    std::map<BinningValue, ProtoLayer::BinningRange> gBinning = {};
    if (not binning.empty()) {
      gBinning = gridBinning(gctx, clusters.second, binning);
    }
    protoLayers.push_back(ProtoLayer(gctx, clusters.second, gBinning));
  }
  return protoLayers;
}

std::vector<Acts::ProtoLayer> Acts::ProtoLayerHelper::protoLayers(
    const GeometryContext& gctx, const std::vector<const Surface*>& surfaces,
    const std::vector<SortingConfig>& sortings,
    const std::vector<BinningConfig>& binning) const {
  ACTS_DEBUG("Received " << surfaces.size() << " surfaces at input.");
  // Finally sorted surfaces
  std::vector<std::vector<const Surface*>> sortSurfaces = {surfaces};
  std::vector<std::map<BinningValue, ProtoLayer::BinningRange>> gridBinnings =
      {};

  for (const auto& sorting : sortings) {
    ACTS_VERBOSE("-> Sorting a set of " << sortSurfaces.size() << " in "
                                        << binningValueNames[sorting.first]);

    std::vector<std::vector<const Surface*>> subSurfaces = {};
    std::vector<std::map<BinningValue, ProtoLayer::BinningRange>> subBinnings =
        {};

    for (const auto& ssurfaces : sortSurfaces) {
      ACTS_VERBOSE("-> Surfaces for this sorting step: " << ssurfaces.size());
      auto pLayers = protoLayers(gctx, ssurfaces, sorting, binning);
      ACTS_VERBOSE("-> Resulted in " << pLayers.size() << " ProtoLayers.");
      for (const auto& pLayer : pLayers) {
        ACTS_VERBOSE("--> ProtoLayer containes " << pLayer.surfaces().size()
                                                 << " surfaces.");
        subSurfaces.push_back(pLayer.surfaces());
        subBinnings.push_back(pLayer.binning);
      }
    }
    sortSurfaces = subSurfaces;
    gridBinnings = subBinnings;
  }
  ACTS_DEBUG("Yielded " << sortSurfaces.size() << " at output.");

  std::vector<Acts::ProtoLayer> finalProtoLayers;

  size_t iss = 0;
  for (auto ssurfaces : sortSurfaces) {
    finalProtoLayers.push_back(
        ProtoLayer(gctx, ssurfaces, gridBinnings[iss++]));
  }

  return finalProtoLayers;
}

std::map<Acts::BinningValue, Acts::ProtoLayer::BinningRange>
Acts::ProtoLayerHelper::gridBinning(
    const GeometryContext& gctx, const std::vector<const Surface*> surfaces,
    const std::vector<BinningConfig>& binning) const {
  /// Initialize
  std::map<BinningValue, std::vector<BinnedRange>> surfaceRanges;
  std::map<BinningValue, BinnedRange> deltaRange;
  for (const auto& bin : binning) {
    ACTS_DEBUG("-> Grid binning in " << binningValueNames[bin.first]
                                     << " with tolerance " << bin.second);
    surfaceRanges[bin.first] = {};
    deltaRange[bin.first] = {std::numeric_limits<double>::max(),
                             -std::numeric_limits<double>::max()};
  }

  for (auto& sf : surfaces) {
    auto sfExtent = sf->polyhedronRepresentation(gctx, 1).extent();
    for (const auto& bin : binning) {
      auto& dtRange = deltaRange[bin.first];
      auto& sfRanges = surfaceRanges[bin.first];
      double bmin = sfExtent.min(bin.first);
      double bmax = sfExtent.max(bin.first);
      ACTS_VERBOSE("--> Bin " << binningValueNames[bin.first] << " checking ["
                              << bmin << ", " << bmax << "].");
      dtRange.first = std::min(bmin, dtRange.first);
      dtRange.second = std::max(bmax, dtRange.second);
      double btol = bin.second;
      if (sfRanges.empty()) {
        sfRanges.push_back({bmin, bmax});
        ACTS_VERBOSE("--> Accepted as new registrant.");
        continue;
      }
      // Check against (all) existing ranges & eventually extend
      bool newRegistrant = true;
      double nrmin = 0.;
      double nrmax = 0.;
      for (auto& sfR : sfRanges) {
        nrmin = sfR.first;
        nrmax = sfR.second;
        if (std::abs(nrmin - bmin) < btol and
            std::abs(sfR.second - nrmax) < btol) {
          newRegistrant = false;
          sfR.first = std::min(nrmin, bmin);
          sfR.second = std::max(nrmax, bmax);
          break;
        }
      }
      if (newRegistrant) {
        sfRanges.push_back({bmin, bmax});
        ACTS_VERBOSE("--> Accepted as new registrant.");
      }
    }
  }

  ACTS_DEBUG("-> Result of the automated binning estimation: ");
  std::map<BinningValue, ProtoLayer::BinningRange> binnings;
  for (const auto& sr : surfaceRanges) {
    ACTS_VERBOSE("--> Binned with " << binningValueNames[sr.first] << " has "
                                    << sr.second.size() << " bins within ["
                                    << deltaRange[sr.first].first << ", "
                                    << deltaRange[sr.first].second << "]");
    binnings[sr.first] =
        ProtoLayer::BinningRange{sr.second.size(), deltaRange[sr.first]};
  }
  return binnings;
}
