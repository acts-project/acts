// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace Acts::Experimental {

/// Maximum number of segments per node
constexpr std::int32_t MAX_SEG_PER_NODE = 1000;
/// Number of segment connections
constexpr std::int32_t N_SEG_CONNS = 6;

class GbtsGeometry;

/// Machine learning lookup table for Gbts seeding
using GbtsMLLookupTable = std::vector<std::array<float, 5>>;

/// GBTs graph node storing space point properties.
class GbtsNode {
 public:
  /// Constructor with layer index
  /// @param l Layer index
  explicit GbtsNode(std::uint16_t l) : m_layer(l) {};
  // readable accessor
  /// Get x coordinate
  /// @return X coordinate
  float x() const noexcept { return m_x; }
  /// Get y coordinate
  /// @return Y coordinate
  float y() const noexcept { return m_y; }
  /// Get phi coordinate
  /// @return Phi coordinate
  float phi() const noexcept { return m_phi; }
  /// Get z coordinate
  /// @return Z coordinate
  float z() const noexcept { return m_z; }
  /// Get radius
  /// @return Radius
  float r() const noexcept { return m_r; }
  /// Get layer index
  /// @return Layer index
  std::uint16_t layer() const noexcept { return m_layer; }
  /// Get pixel cluster width
  /// @return Pixel cluster width
  float pixelClusterWidth() const noexcept { return m_pcw; }
  /// Get local position y
  /// @return Local y position
  float localPositionY() const noexcept { return m_locPosY; }
  /// Get space point index
  /// @return Space point index
  std::uint32_t sp_idx() const noexcept { return m_idx; }

  // writeable accessor
  /// Get mutable x coordinate
  /// @return Reference to x coordinate
  float& x() noexcept { return m_x; }
  /// Get mutable y coordinate
  /// @return Reference to y coordinate
  float& y() noexcept { return m_y; }
  /// Get mutable phi coordinate
  /// @return Reference to phi coordinate
  float& phi() noexcept { return m_phi; }
  /// Get mutable z coordinate
  /// @return Reference to z coordinate
  float& z() noexcept { return m_z; }
  /// Get mutable radius
  /// @return Reference to radius
  float& r() noexcept { return m_r; }
  /// Get mutable pixel cluster width
  /// @return Reference to pixel cluster width
  float& pixelClusterWidth() noexcept { return m_pcw; }
  /// Get mutable local position y
  /// @return Reference to local y position
  float& localPositionY() noexcept { return m_locPosY; }
  /// Get mutable space point index
  /// @return Reference to space point index
  std::uint32_t& sp_idx() noexcept { return m_idx; }

 private:
  float m_x{};
  float m_y{};
  float m_z{};
  float m_r{};
  float m_phi{};
  std::uint16_t m_layer{10000};
  std::uint32_t m_idx{std::numeric_limits<std::uint32_t>::max()};
  float m_pcw{};
  float m_locPosY{};
};

/// Eta-bin container for GBTs nodes and edge data.
class GbtsEtaBin {
 public:
  /// Comparator sorting nodes by increasing phi.
  struct CompareNodesByPhi {
    /// Compare two nodes by phi
    /// @param n1 First node
    /// @param n2 Second node
    /// @return True if first node has smaller phi
    bool operator()(const GbtsNode* n1, const GbtsNode* n2) {
      return n1->phi() < n2->phi();
    }
  };

  GbtsEtaBin();

  /// Sort nodes by phi
  void sortByPhi();
  /// Initialize node attributes
  void initializeNodes();
  /// Check if bin is empty
  /// @return True if bin has no nodes
  bool empty() const { return m_vn.empty(); }

  /// Generate phi indexing
  /// @param dphi Phi bin width
  void generatePhiIndexing(float dphi);

  /// Get minimum bin radius
  /// @return Minimum radius in this bin
  float getMinBinRadius() const { return m_minRadius; }

  /// Get maximum bin radius
  /// @return Maximum radius in this bin
  float getMaxBinRadius() const { return m_maxRadius; }

  /// nodes of the graph
  std::vector<const GbtsNode*> m_vn;
  /// Phi-indexed nodes
  std::vector<std::pair<float, std::uint32_t>> m_vPhiNodes;
  /// vectors of incoming edges, stores indices of edges in the edge vector
  std::vector<std::vector<std::uint32_t>> m_in;
  /// node attributes: m_minCutOnTau, m_maxCutOnTau, m_phi, m_r, m_z;
  std::vector<std::array<float, 5>> m_params;
  /// Minimum radius in bin
  float m_minRadius{};
  /// Maximum radius in bin
  float m_maxRadius{};

  /// Layer key for this bin
  std::uint32_t m_layerKey{0};
};

/// Storage container for GBTs nodes and edges.
class GbtsDataStorage {
 public:
  /// Constructor
  /// @param geometry Shared pointer to GBTs geometry
  /// @param config Configuration for seed finder
  /// @param mlLUT Machine learning lookup table
  explicit GbtsDataStorage(std::shared_ptr<const GbtsGeometry> geometry,
                           const SeedFinderGbtsConfig& config,
                           GbtsMLLookupTable mlLUT);

  /// Load pixel graph nodes
  /// @param layerIndex Layer index for the nodes
  /// @param coll Collection of nodes to load
  /// @param useML Use machine learning features
  /// @return Number of nodes loaded
  std::int32_t loadPixelGraphNodes(std::int16_t layerIndex,
                                   const std::span<const GbtsNode> coll,
                                   bool useML);
  /// Load strip graph nodes
  /// @param layerIndex Layer index for the nodes
  /// @param coll Collection of nodes to load
  /// @return Number of nodes loaded
  std::int32_t loadStripGraphNodes(std::int16_t layerIndex,
                                   const std::span<const GbtsNode> coll);

  /// Get the total number of nodes
  /// @return Total number of nodes
  std::uint32_t numberOfNodes() const;
  /// Sort nodes by phi
  void sortByPhi();
  /// Initialize node attributes
  /// @param useML Use machine learning features
  void initializeNodes(bool useML);
  /// Generate phi indexing
  /// @param dphi Phi bin width
  void generatePhiIndexing(float dphi);

  /// Get eta bin by index
  /// @param idx Eta bin index
  /// @return Reference to the eta bin
  GbtsEtaBin& getEtaBin(std::int32_t idx) {
    if (idx >= static_cast<std::int32_t>(m_etaBins.size())) {
      idx = idx - 1;
    }
    return m_etaBins.at(idx);
  }

 protected:
  /// GBTs geometry
  std::shared_ptr<const GbtsGeometry> m_geo;

  /// Configuration for seed finder
  SeedFinderGbtsConfig m_config;

  /// Machine learning lookup table
  GbtsMLLookupTable m_mlLUT;

  /// Eta bins for node storage
  std::vector<GbtsEtaBin> m_etaBins;
};

/// Edge between two GBTs nodes with fit parameters.
class GbtsEdge {
 public:
  /// Comparator sorting edges by descending level.
  struct CompareLevel {
   public:
    /// Compare two edges by level
    /// @param pE1 First edge
    /// @param pE2 Second edge
    /// @return True if first edge has higher level than second
    bool operator()(const GbtsEdge* pE1, const GbtsEdge* pE2) {
      return pE1->m_level > pE2->m_level;
    }
  };

  GbtsEdge() = default;

  /// Constructor
  /// @param n1 First node
  /// @param n2 Second node
  /// @param p1 First fit parameter
  /// @param p2 Second fit parameter
  /// @param p3 Third fit parameter
  GbtsEdge(const GbtsNode* n1, const GbtsNode* n2, float p1, float p2, float p3)
      : m_n1(n1), m_n2(n2), m_level(1), m_next(1) {
    m_p[0] = p1;
    m_p[1] = p2;
    m_p[2] = p3;
  }

  /// First node of the edge
  const GbtsNode* m_n1{nullptr};
  /// Second node of the edge
  const GbtsNode* m_n2{nullptr};

  /// Level in the graph hierarchy
  std::int8_t m_level{-1};
  /// Index of next edge
  std::int8_t m_next{-1};

  /// Number of neighbor edges
  std::uint8_t m_nNei{0};
  /// Fit parameters
  std::array<float, 3> m_p{};

  /// Global indices of the connected edges
  std::array<std::uint32_t, N_SEG_CONNS> m_vNei{};
};

}  // namespace Acts::Experimental
