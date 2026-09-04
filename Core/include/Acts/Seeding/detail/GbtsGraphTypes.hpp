// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/StripSpacePointCalibrationDetails.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

namespace Acts::Experimental::detail {

/// Accepted |cot(theta)| range for one cluster width bin. A cluster near the
/// module edge may be shortened, hence the second pair.
struct GbtsTauBounds final {
  /// Minimum accepted |cot(theta)|.
  float minTau{};
  /// Maximum accepted |cot(theta)|. Negative means undertrained: do not cut.
  float maxTau{};
  /// Minimum accepted |cot(theta)| near the module edge.
  float minTauNearEdge{};
  /// Maximum accepted |cot(theta)| near the module edge. Negative as above.
  float maxTauNearEdge{};
};

/// Tau bounds per cluster width, one entry per 0.05 mm, indexed not searched.
using GbtsTauLookupTable = std::vector<GbtsTauBounds>;

/// Maximum number of neighbouring edges recorded per graph edge
static constexpr std::uint32_t kGbtsMaxEdgeNeighbours = 6;

/// Bins of the per-node z0 histogram, which is kept as a bit mask in
/// GbtsNodeEdgeInfo::isConnected and so may not exceed its width.
static constexpr std::int32_t kGbtsZ0HistogramBins = 16;
static_assert(kGbtsZ0HistogramBins <=
              std::numeric_limits<std::uint16_t>::digits);

//! [gbts node params]
/// Per-node parameters used while building the graph.
///
/// All five values are read together in the innermost doublet loop, so they sit
/// in one record rather than in separate arrays.
struct GbtsNodeParams final {
  /// Minimum accepted |cot(theta)|. The infinite defaults mean "do not cut on
  /// tau"; only the tau lookup table narrows them.
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
//! [gbts node params]

//! [gbts node edge info]
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
//! [gbts node edge info]

/// One outer eta bin reachable from an inner one, with the cut the layer pair
/// they belong to carries.
struct GbtsBinLink final {
  /// The outer eta bin.
  std::uint32_t bin{};
  /// GbtsLayerConnection::tauRatioCut of the pair this link came from.
  float tauRatioCut{};
};

//! [gbts eta bin info]
/// Constant per-eta-bin data.
struct GbtsEtaBinInfo final {
  /// Range of node indices belonging to this bin.
  SpacePointIndexRange nodes{0, 0};

  /// (phi, node index) pairs covering this bin, with duplicated entries shifted
  /// by +-2pi at the wrap-around so the phi sliding window never has to wrap.
  std::vector<std::pair<float, SpacePointIndex>> phiNodes;

  float minRadius{};
  float maxRadius{};
  GbtsExperimentLayerId layerId{0};

  /// Type of the layer this bin belongs to.
  GbtsLayerType type{};
  /// Technology of the layer this bin belongs to.
  GbtsLayerTechnology technology{};

  /// Check if bin is empty
  /// @return True if bin has no nodes
  bool empty() const { return nodes.first == nodes.second; }
};
//! [gbts eta bin info]

class GbtsNodeProxy;

/// What a node's entry in `GbtsNodeView::stripIndex` holds when it carries no
/// stereo pair, i.e. when it is a pixel node.
constexpr std::uint32_t kNoStrip = std::numeric_limits<std::uint32_t>::max();

/// Read-only view of the node attributes needed outside the graph builder.
///
/// Positions are packed as (x, y, z, r) because every consumer reads them
/// together. Index it to get a GbtsNodeProxy with named accessors.
struct GbtsNodeView final {
  /// Packed (x, y, z, r) per node.
  std::span<const std::array<float, 4>> positions;
  /// Dense layer index per node, narrowed to 16 bits.
  std::span<const std::uint16_t> layers;
  /// Stereo pairs of the strip nodes, reached through `stripIndex`.
  std::span<const OuterStripSpacePointCalibrationDetailsDerived> strips;
  /// Index into `strips` per node, `kNoStrip` where a node carries none.
  /// Empty when nothing does.
  std::span<const std::uint32_t> stripIndex;

  /// Handle on a single node.
  /// @param index Index of the node
  /// @return Proxy for the node
  GbtsNodeProxy operator[](SpacePointIndex index) const;

  /// The stereo pair a node was formed from.
  /// @param index Index of the node
  /// @return The pair, or null for a pixel node
  const OuterStripSpacePointCalibrationDetailsDerived* strip(
      SpacePointIndex index) const {
    if (stripIndex.empty() || stripIndex[index] == kNoStrip) {
      return nullptr;
    }
    return &strips[stripIndex[index]];
  }
};

/// Read-only handle on a single graph node.
///
/// Cheap to copy and meant to be created at the point of use. The view it
/// refers to has to outlive it.
class GbtsNodeProxy final {
 public:
  /// @param view View of the node positions and layers
  /// @param index Index of the node
  GbtsNodeProxy(const GbtsNodeView& view, SpacePointIndex index)
      : m_view(&view), m_index(index) {}

  /// @return Index of the node
  SpacePointIndex index() const { return m_index; }
  /// @return Global x coordinate
  float x() const { return position()[0]; }
  /// @return Global y coordinate
  float y() const { return position()[1]; }
  /// @return Global z coordinate
  float z() const { return position()[2]; }
  /// @return Transverse distance from the beamline
  float r() const { return position()[3]; }
  /// @return Dense layer index
  GbtsLayerIndex layer() const { return m_view->layers[m_index]; }

 private:
  const std::array<float, 4>& position() const {
    return m_view->positions[m_index];
  }

  const GbtsNodeView* m_view{};
  SpacePointIndex m_index{kSpacePointIndexInvalid};
};

inline GbtsNodeProxy GbtsNodeView::operator[](SpacePointIndex index) const {
  return GbtsNodeProxy(*this, index);
}

//! [gbts edge]
/// Edge between two GBTS nodes with fit parameters.
struct GbtsEdge final {
  GbtsEdge() = default;

  /// Constructor
  /// @param n1_ Inner node index
  /// @param n2_ Outer node index
  /// @param tauRatioCut_ Cut the layer pair of this doublet carries
  /// @param n2PixelBarrel_ Whether the outer node is on a pixel barrel layer
  /// @param p1_ First fit parameter
  /// @param p2_ Second fit parameter
  /// @param p3_ Third fit parameter
  GbtsEdge(SpacePointIndex n1_, SpacePointIndex n2_, float tauRatioCut_,
           bool n2PixelBarrel_, float p1_, float p2_, float p3_)
      : n1{n1_},
        n2{n2_},
        level{1},
        next{1},
        n2PixelBarrel{n2PixelBarrel_},
        p{p1_, p2_, p3_},
        tauRatioCut{tauRatioCut_} {}

  /// Inner node of the edge
  SpacePointIndex n1{kSpacePointIndexInvalid};
  /// Outer node of the edge
  SpacePointIndex n2{kSpacePointIndexInvalid};

  std::int8_t level{-1};
  std::int8_t next{-1};

  std::uint8_t nNei{0};

  /// Whether the outer node is on a pixel barrel layer, the only thing the
  /// innermost neighbour loop asks about it. Cached so that loop does not have
  /// to chase the node's bin, and in what was padding so the edge does not
  /// grow.
  bool n2PixelBarrel{};

  std::array<float, 3> p{};

  /// GbtsLayerConnection::tauRatioCut of this doublet's layer pair. Cached
  /// next to the fit parameters so the innermost neighbour loop does not have
  /// to chase the pair.
  float tauRatioCut{};

  std::array<std::uint32_t, kGbtsMaxEdgeNeighbours> vNei{};
};
//! [gbts edge]

}  // namespace Acts::Experimental::detail
