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

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace Acts::Experimental {

/// Maximum number of neighbouring edges recorded per graph edge
static constexpr std::uint32_t kGbtsMaxEdgeNeighbours = 6;

/// Number of values per row of the machine learning lookup table: the cluster
/// width the row is keyed on, followed by two pairs of tau bounds.
static constexpr std::size_t kGbtsMlLookupTableColumns = 5;

class GbtsGeometry;

/// Machine learning lookup table for Gbts seeding
using GbtsMlLookupTable =
    std::vector<std::array<float, kGbtsMlLookupTableColumns>>;

/// Index of a graph node inside GbtsNodeStorage. Nodes are stored ordered by
/// (eta bin, phi), so all nodes of an eta bin form a contiguous range.
using GbtsNodeIndex = SpacePointIndex;

/// Sentinel for an unset graph node index.
static constexpr GbtsNodeIndex kGbtsNodeIndexInvalid = kSpacePointIndexInvalid;

/// Per-node parameters used while building the graph.
///
/// All five values are read together in the innermost doublet loop, so they sit
/// in one record rather than in separate arrays.
struct GbtsNodeParams final {
  /// Minimum accepted |cot(theta)|. The infinite defaults mean "do not cut on
  /// tau"; only the machine learning lookup table narrows them.
  float minTau{-std::numeric_limits<float>::infinity()};
  /// Maximum accepted |cot(theta)|. See @ref minTau.
  float maxTau{std::numeric_limits<float>::infinity()};
  /// Azimuthal angle in the xy plane.
  float phi{};
  /// Transverse distance from the beamline.
  float r{};
  /// Global z coordinate.
  float z{};
};

/// Per-node graph bookkeeping, written while the graph is built.
///
/// The three fields are read together in the innermost doublet loop, so they
/// sit in one record rather than in separate arrays.
struct GbtsNodeEdgeInfo final {
  /// Index of the first incoming graph edge attached to the node.
  std::uint32_t firstEdge{};
  /// Total number of incoming graph edges attached to the node.
  std::uint16_t numEdges{};
  /// z0 histogram bitmask. Non-zero marks the node's outer neighbourhood as
  /// connected to the previously built graph.
  std::uint16_t isConnected{};
};

/// Constant per-eta-bin data.
struct GbtsEtaBinInfo final {
  /// Range of node indices belonging to this bin.
  SpacePointIndexRange nodes{0, 0};

  /// (phi, node index) pairs covering this bin, with duplicated entries shifted
  /// by +-2pi at the wrap-around so the phi sliding window never has to wrap.
  std::vector<std::pair<float, GbtsNodeIndex>> phiNodes;

  /// Minimum radius in bin
  float minRadius{};
  /// Maximum radius in bin
  float maxRadius{};
  /// GBTS layer ID for this bin
  std::uint32_t layerId{0};

  /// Check if bin is empty
  /// @return True if bin has no nodes
  bool empty() const { return nodes.first == nodes.second; }
};

class GbtsNodeProxy;

/// Read-only view of the node attributes needed outside the graph builder.
///
/// Positions are packed as (x, y, z, r) because every consumer reads them
/// together. Index it to get a GbtsNodeProxy with named accessors.
struct GbtsNodeView final {
  /// Packed (x, y, z, r) per node.
  std::span<const std::array<float, 4>> positions;
  /// Dense layer index per node.
  std::span<const std::uint16_t> layers;

  /// Handle on a single node.
  /// @param index Index of the node
  /// @return Proxy for the node
  GbtsNodeProxy operator[](GbtsNodeIndex index) const;
};

/// Read-only handle on a single graph node.
///
/// Cheap to copy and meant to be created at the point of use. The view it
/// refers to has to outlive it.
class GbtsNodeProxy final {
 public:
  /// @param view View of the node positions and layers
  /// @param index Index of the node
  GbtsNodeProxy(const GbtsNodeView& view, GbtsNodeIndex index)
      : m_view(&view), m_index(index) {}

  /// @return Index of the node
  GbtsNodeIndex index() const { return m_index; }
  /// @return Global x coordinate
  float x() const { return position()[0]; }
  /// @return Global y coordinate
  float y() const { return position()[1]; }
  /// @return Global z coordinate
  float z() const { return position()[2]; }
  /// @return Transverse distance from the beamline
  float r() const { return position()[3]; }
  /// @return Dense layer index
  std::uint16_t layer() const { return m_view->layers[m_index]; }

 private:
  const std::array<float, 4>& position() const {
    return m_view->positions[m_index];
  }

  const GbtsNodeView* m_view{};
  GbtsNodeIndex m_index{kGbtsNodeIndexInvalid};
};

inline GbtsNodeProxy GbtsNodeView::operator[](GbtsNodeIndex index) const {
  return GbtsNodeProxy(*this, index);
}

/// Edge between two GBTS nodes with fit parameters.
struct GbtsEdge final {
  GbtsEdge() = default;

  /// Constructor
  /// @param n1_ Inner node index
  /// @param n2_ Outer node index
  /// @param n2LayerId_ GBTS layer ID of the outer node
  /// @param p1_ First fit parameter
  /// @param p2_ Second fit parameter
  /// @param p3_ Third fit parameter
  GbtsEdge(GbtsNodeIndex n1_, GbtsNodeIndex n2_, std::uint32_t n2LayerId_,
           float p1_, float p2_, float p3_)
      : n1{n1_},
        n2{n2_},
        level{1},
        next{1},
        p{p1_, p2_, p3_},
        n2LayerId{n2LayerId_} {}

  /// Inner node of the edge
  GbtsNodeIndex n1{kGbtsNodeIndexInvalid};
  /// Outer node of the edge
  GbtsNodeIndex n2{kGbtsNodeIndexInvalid};

  /// Level in the graph hierarchy
  std::int8_t level{-1};
  /// Index of next edge
  std::int8_t next{-1};

  /// Number of neighbor edges
  std::uint8_t nNei{0};
  /// Fit parameters
  std::array<float, 3> p{};

  /// GBTS layer ID of the outer node. Cached next to the fit parameters so the
  /// innermost neighbour loop does not have to chase the node.
  std::uint32_t n2LayerId{0};

  /// Global indices of the connected edges
  std::array<std::uint32_t, kGbtsMaxEdgeNeighbours> vNei{};
};

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
    /// Use machine learning features (cluster width based tau cuts).
    bool useMl = false;
    /// Maximum endcap cluster width, applied to pixel endcap nodes when
    /// machine learning features are enabled.
    float maxEndcapClusterWidth = 0.35f;
    /// Half-length in local y of a pixel module, against which the distance of
    /// a cluster to the module edge is measured.
    float moduleHalfLengthY = 10.f;
    /// Distance to the module edge below which a cluster may be shortened. The
    /// machine learning lookup table carries a separate pair of tau bounds for
    /// that case.
    float moduleEdgeTolerance = 0.3f;
    /// Width of the phi slice used to build the phi indexing.
    float phiSliceWidth = 0.f;
  };

  /// @param config Node loading configuration
  /// @param geometry Shared pointer to GBTS geometry
  /// @param mlLut Machine learning lookup table
  GbtsNodeStorage(Config config, std::shared_ptr<const GbtsGeometry> geometry,
                  GbtsMlLookupTable mlLut);

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

  /// Get eta bin info by index
  /// @param idx Eta bin index
  /// @return Reference to the eta bin info
  const GbtsEtaBinInfo& etaBin(std::uint32_t idx) const {
    return m_etaBins.at(idx < m_etaBins.size() ? idx : idx - 1);
  }

  /// Read-only view of the node positions and layers
  /// @return Node view
  GbtsNodeView nodeView() const {
    return GbtsNodeView{m_nodes.xyzrColumn().data(), m_layers};
  }

  /// Per-node graph parameters, indexed by node index
  /// @return Span over the node parameters
  std::span<const GbtsNodeParams> nodeParams() const {
    return m_paramsColumn->data();
  }

  /// Per-node graph bookkeeping, indexed by node index
  /// @return Mutable span over the node edge info
  std::span<GbtsNodeEdgeInfo> nodeEdgeInfo() {
    return m_edgeInfoColumn->data();
  }

  /// Per-node graph bookkeeping, indexed by node index
  /// @return Span over the node edge info
  std::span<const GbtsNodeEdgeInfo> nodeEdgeInfo() const {
    return m_edgeInfoColumn->data();
  }

  /// Map a node index back to the index the caller used when inserting it.
  /// @param node Node index
  /// @return The caller's space point index
  SpacePointIndex spacePointIndex(GbtsNodeIndex node) const {
    return m_nodes.copiedFromIndexColumn()[node];
  }

 private:
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

  /// Narrow a node's tau window using the machine learning lookup table.
  /// @param staged The staged node
  /// @param params The node parameters to narrow
  void applyMlTauCuts(const StagedNode& staged, GbtsNodeParams& params) const;

  /// Build the wrap-around aware phi indexing for every bin.
  /// @param dphi Width of the phi margin duplicated at the wrap-around
  void generatePhiIndexing(float dphi);

  Config m_cfg;

  std::shared_ptr<const GbtsGeometry> m_geometry;

  GbtsMlLookupTable m_mlLut;

  /// Nodes ordered by (eta bin, phi). Carries the caller's index and the packed
  /// (x, y, z, r) position, plus the derived data as dynamic columns.
  SpacePointContainer m_nodes;

  std::optional<MutableSpacePointColumnProxy<GbtsNodeParams>> m_paramsColumn;
  std::optional<MutableSpacePointColumnProxy<GbtsNodeEdgeInfo>>
      m_edgeInfoColumn;

  /// Dense layer index per node, in node order.
  std::vector<std::uint16_t> m_layers;

  std::vector<GbtsEtaBinInfo> m_etaBins;

  /// Nodes as inserted, before sorting.
  std::vector<StagedNode> m_staged;
  /// Staged node indices per eta bin.
  std::vector<std::vector<std::uint32_t>> m_stagedPerBin;
};

}  // namespace Acts::Experimental
