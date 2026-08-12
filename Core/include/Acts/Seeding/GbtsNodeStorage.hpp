// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointColumnProxy.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/GbtsTypes.hpp"
#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace Acts::Experimental {

class GraphBasedTrackSeeder;
class GbtsGeometry;

/// Storage for the GBTS graph nodes.
///
/// Nodes go in one at a time through `insert`, which takes plain scalars so
/// that a caller can fill the storage from its own space point EDM. `finalize`
/// then orders the nodes by (eta bin, phi) into a space point container, with
/// the derived per-node data in dynamic columns on that container.
class GbtsNodeStorage final {
 public:
  /// Configuration for node loading.
  struct Config {
    /// Per-layer flag marking pixel layers. Indexed by dense layer index.
    std::vector<bool> isPixelLayer;
    /// Cut on the pixel cluster width: wide endcap clusters and the tau
    /// lookup table.
    bool useClusterWidthCuts = false;
    /// Maximum endcap cluster width, applied to pixel endcap nodes when
    /// the cluster width cuts are enabled.
    float maxEndcapClusterWidth = 0.35f;
    /// Half-length in local y of a pixel module, against which the distance of
    /// a cluster to the module edge is measured.
    float moduleHalfLengthY = 10.f;
    /// Distance to the module edge below which a cluster may be shortened,
    /// which the tau lookup table covers with its own bounds.
    float moduleEdgeTolerance = 0.3f;
    /// Width of the phi slice used to build the phi indexing.
    float phiSliceWidth = 0.f;
  };

  /// @param config Node loading configuration
  /// @param geometry Shared pointer to GBTS geometry
  /// @param tauLut Per-cluster-width tau bounds
  GbtsNodeStorage(Config config, std::shared_ptr<const GbtsGeometry> geometry,
                  GbtsTauLookupTable tauLut);

  /// Insert a space point, deriving r and phi from the global position.
  /// @param index Index of the space point in the caller's own collection
  /// @param x Global x coordinate
  /// @param y Global y coordinate
  /// @param z Global z coordinate
  /// @param layerIndex Dense GBTS layer index
  /// @param clusterWidth Pixel cluster width
  /// @param localPositionY Local y cluster position
  /// @return The eta bin the node was placed in, or nullopt if it was rejected
  std::optional<std::uint32_t> insert(SpacePointIndex index, float x, float y,
                                      float z, std::uint32_t layerIndex,
                                      float clusterWidth = 0.f,
                                      float localPositionY = 0.f);

  /// Insert a space point for callers that already have r and phi.
  /// @param index Index of the space point in the caller's own collection
  /// @param x Global x coordinate
  /// @param y Global y coordinate
  /// @param z Global z coordinate
  /// @param r Transverse distance from the beamline
  /// @param phi Azimuthal angle in the xy plane
  /// @param layerIndex Dense GBTS layer index
  /// @param clusterWidth Pixel cluster width
  /// @param localPositionY Local y cluster position
  /// @return The eta bin the node was placed in, or nullopt if it was rejected
  std::optional<std::uint32_t> insert(SpacePointIndex index, float x, float y,
                                      float z, float r, float phi,
                                      std::uint32_t layerIndex,
                                      float clusterWidth = 0.f,
                                      float localPositionY = 0.f);

  /// Insert a space point from an ACTS space point container.
  /// @param sp The space point to insert
  /// @param layerColumn Column holding the dense GBTS layer index
  /// @param clusterWidthColumn Column holding the pixel cluster width
  /// @param localPositionYColumn Column holding the local y cluster position
  /// @return The eta bin the node was placed in, or nullopt if it was rejected
  std::optional<std::uint32_t> insert(
      const ConstSpacePointProxy& sp,
      const ConstSpacePointColumnProxy<std::uint32_t>& layerColumn,
      const ConstSpacePointColumnProxy<float>& clusterWidthColumn,
      const ConstSpacePointColumnProxy<float>& localPositionYColumn) {
    return insert(sp.index(), sp.x(), sp.y(), sp.z(), sp.r(), sp.phi(),
                  sp.extra(layerColumn), sp.extra(clusterWidthColumn),
                  sp.extra(localPositionYColumn));
  }

  /// Insert every space point of a container.
  /// @param spacePoints The space points to insert
  /// @param layerColumn Column holding the dense GBTS layer index
  /// @param clusterWidthColumn Column holding the pixel cluster width
  /// @param localPositionYColumn Column holding the local y cluster position
  void extend(const SpacePointContainer& spacePoints,
              const ConstSpacePointColumnProxy<std::uint32_t>& layerColumn,
              const ConstSpacePointColumnProxy<float>& clusterWidthColumn,
              const ConstSpacePointColumnProxy<float>& localPositionYColumn);

  /// Sort the nodes by (eta bin, phi) and build the derived per-node data.
  /// Must be called once after all inserts and before the storage is read.
  void finalize();

  /// Get the total number of nodes
  /// @return Total number of nodes
  std::uint32_t numberOfNodes() const { return m_nodes.size(); }

  /// Map a node index back to the index the caller used when inserting it.
  /// @param node Node index
  /// @return The caller's space point index
  SpacePointIndex spacePointIndex(GbtsNodeIndex node) const {
    return m_nodes.copiedFromIndexColumn()[node];
  }

 private:
  // The graph representation is an implementation detail shared only with the
  // seeder that walks it.
  friend class GraphBasedTrackSeeder;

  /// Get eta bin info by index
  /// @param idx Eta bin index
  /// @return Reference to the eta bin info
  const detail::GbtsEtaBinInfo& etaBin(std::uint32_t idx) const {
    return m_etaBins.at(idx < m_etaBins.size() ? idx : idx - 1);
  }

  /// Read-only view of the node positions and layers
  /// @return Node view
  detail::GbtsNodeView nodeView() const {
    return detail::GbtsNodeView{m_nodes.xyzrColumn().data(), m_layers};
  }

  /// Per-node graph parameters, indexed by node index
  /// @return Span over the node parameters
  std::span<const detail::GbtsNodeParams> nodeParams() const {
    return m_paramsColumn->data();
  }

  /// Per-node graph bookkeeping, indexed by node index
  /// @return Mutable span over the node edge info
  std::span<detail::GbtsNodeEdgeInfo> nodeEdgeInfo() {
    return m_edgeInfoColumn->data();
  }

  /// Per-node graph bookkeeping, indexed by node index
  /// @return Span over the node edge info
  std::span<const detail::GbtsNodeEdgeInfo> nodeEdgeInfo() const {
    return m_edgeInfoColumn->data();
  }

  /// A node as recorded by `insert`, before sorting.
  struct StagedNode {
    SpacePointIndex spacePointIndex{};
    float x{};
    float y{};
    float z{};
    float r{};
    float phi{};
    float clusterWidth{};
    float localPositionY{};
    std::uint16_t layer{};
  };

  /// Sort a single bin's staged nodes by phi.
  /// @param staged Staged node indices of the bin
  /// @return The staged indices in phi order
  std::vector<std::uint32_t> sortBinByPhi(
      const std::vector<std::uint32_t>& staged) const;

  /// Narrow a node's tau window using the tau lookup table.
  /// @param staged The staged node
  /// @param params The node parameters to narrow
  void applyTauCuts(const StagedNode& staged,
                    detail::GbtsNodeParams& params) const;

  /// Build the wrap-around aware phi indexing for every bin.
  /// @param dphi Width of the phi margin duplicated at the wrap-around
  void generatePhiIndexing(float dphi);

  Config m_cfg;

  std::shared_ptr<const GbtsGeometry> m_geometry;

  GbtsTauLookupTable m_tauLut;

  /// Nodes ordered by (eta bin, phi). Carries the caller's index and the packed
  /// (x, y, z, r) position, plus the derived data as dynamic columns.
  SpacePointContainer m_nodes;

  std::optional<MutableSpacePointColumnProxy<detail::GbtsNodeParams>>
      m_paramsColumn;
  std::optional<MutableSpacePointColumnProxy<detail::GbtsNodeEdgeInfo>>
      m_edgeInfoColumn;

  /// Dense layer index per node, in node order.
  std::vector<std::uint16_t> m_layers;

  std::vector<detail::GbtsEtaBinInfo> m_etaBins;

  /// Nodes as inserted, before sorting.
  std::vector<StagedNode> m_staged;
  /// Staged node indices per eta bin.
  std::vector<std::vector<std::uint32_t>> m_stagedPerBin;
};

}  // namespace Acts::Experimental
