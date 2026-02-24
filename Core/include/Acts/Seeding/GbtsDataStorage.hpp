// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsConfig.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <vector>

namespace Acts::Experimental {

/// Number of segment connections
constexpr std::uint32_t gbtsNumSegConns = 6;

class GbtsGeometry;

/// Machine learning lookup table for Gbts seeding
using GbtsMlLookupTable = std::vector<std::array<float, 5>>;

/// GBTs graph node storing space point properties.
struct GbtsNode final {
 public:
  /// Constructor with layer index
  /// @param layer_ Layer index
  explicit GbtsNode(std::uint16_t layer_) : layer(layer_) {};

  /// Global x coordinate of the node
  float x{};
  /// Global y coordinate of the node
  float y{};
  /// Global z coordinate of the node
  float z{};
  /// Transverse distance from the beamline
  float r{};
  /// Azimuthal angle in the xy plane
  float phi{};
  /// Layer index
  std::uint16_t layer{};
  /// Index of the node in the original collection
  std::uint32_t idx{std::numeric_limits<std::uint32_t>::max()};
  /// Pixel cluster width
  float pcw{};
  /// Local y cluster position
  float locPosY{};
};

/// Eta-bin container for GBTs nodes and edge data.
struct GbtsEtaBin final {
  GbtsEtaBin();

  /// Sort nodes by phi
  void sortByPhi();
  /// Initialize node attributes
  void initializeNodes();
  /// Check if bin is empty
  /// @return True if bin has no nodes
  bool empty() const { return vn.empty(); }

  /// Generate phi indexing
  /// @param dphi Phi bin width
  void generatePhiIndexing(float dphi);

  /// nodes of the graph
  std::vector<const GbtsNode*> vn;
  /// Phi-indexed nodes
  std::vector<std::pair<float, std::uint32_t>> vPhiNodes;
  /// node attributes: minCutOnTau, maxCutOnTau, phi, r, z;
  std::vector<std::array<float, 5>> params;
  /// the index of the first incoming graph edge attached to the node
  std::vector<std::uint32_t> vFirstEdge;
  /// the total number of incoming graph edges attached to this node
  std::vector<std::uint16_t> vNumEdges;
  /// flag to indicate the node's outer neighbourhood isolation from previously
  /// built graph
  std::vector<std::uint16_t> vIsConnected;

  /// Minimum radius in bin
  float minRadius{};
  /// Maximum radius in bin
  float maxRadius{};

  /// Layer key for this bin
  std::uint32_t layerKey{0};
};

/// Storage container for GBTs nodes and edges.
class GbtsDataStorage final {
 public:
  /// Constructor
  /// @param config Configuration for seed finder
  /// @param geometry Shared pointer to GBTs geometry
  /// @param mlLut Machine learning lookup table
  explicit GbtsDataStorage(const GbtsConfig& config,
                           std::shared_ptr<const GbtsGeometry> geometry,
                           GbtsMlLookupTable mlLut);

  /// Load pixel graph nodes
  /// @param layerIndex Layer index for the nodes
  /// @param coll Collection of nodes to load
  /// @param useMl Use machine learning features
  /// @return Number of nodes loaded
  std::uint32_t loadPixelGraphNodes(std::uint16_t layerIndex,
                                    const std::span<const GbtsNode> coll,
                                    bool useMl);
  /// Load strip graph nodes
  /// @param layerIndex Layer index for the nodes
  /// @param coll Collection of nodes to load
  /// @return Number of nodes loaded
  std::uint32_t loadStripGraphNodes(std::uint16_t layerIndex,
                                    const std::span<const GbtsNode> coll);

  /// Get the total number of nodes
  /// @return Total number of nodes
  std::uint32_t numberOfNodes() const;
  /// Sort nodes by phi
  void sortByPhi();
  /// Initialize node attributes
  /// @param useMl Use machine learning features
  void initializeNodes(bool useMl);
  /// Generate phi indexing
  /// @param dphi Phi bin width
  void generatePhiIndexing(float dphi);

  /// Get eta bin by index
  /// @param idx Eta bin index
  /// @return Reference to the eta bin
  GbtsEtaBin& getEtaBin(std::uint32_t idx) {
    if (idx >= m_etaBins.size()) {
      idx = idx - 1;
    }
    return m_etaBins.at(idx);
  }

 private:
  /// GBTs geometry
  std::shared_ptr<const GbtsGeometry> m_geo;

  /// Configuration for seed finder
  GbtsConfig m_cfg{};

  /// Machine learning lookup table
  GbtsMlLookupTable m_mlLut;

  /// Eta bins for node storage
  std::vector<GbtsEtaBin> m_etaBins;
};

/// Edge between two GBTs nodes with fit parameters.
struct GbtsEdge final {
  GbtsEdge() = default;

  /// Constructor
  /// @param n1_ First node
  /// @param n2_ Second node
  /// @param p1_ First fit parameter
  /// @param p2_ Second fit parameter
  /// @param p3_ Third fit parameter
  GbtsEdge(const GbtsNode* n1_, const GbtsNode* n2_, float p1_, float p2_,
           float p3_)
      : n1{n1_}, n2{n2_}, level{1}, next{1}, p{p1_, p2_, p3_} {}

  /// First node of the edge
  const GbtsNode* n1{nullptr};
  /// Second node of the edge
  const GbtsNode* n2{nullptr};

  /// Level in the graph hierarchy
  std::int8_t level{-1};
  /// Index of next edge
  std::int8_t next{-1};

  /// Number of neighbor edges
  std::uint8_t nNei{0};
  /// Fit parameters
  std::array<float, 3> p{};

  /// Global indices of the connected edges
  std::array<std::uint32_t, gbtsNumSegConns> vNei{};
};

}  // namespace Acts::Experimental
