// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>
namespace Acts::Experimental {

constexpr int MAX_SEG_PER_NODE = 1000;
constexpr int N_SEG_CONNS = 6;

class GbtsGeometry;

using GbtsMLLookupTable = std::vector<std::array<float, 5>>;

class GbtsNode {
 public:
  explicit GbtsNode(std::uint16_t l) : m_layer(l) {};
  // readable accessor
  float x() const noexcept { return m_x; }
  float y() const noexcept { return m_y; }
  float phi() const noexcept { return m_phi; }
  float z() const noexcept { return m_z; }
  float r() const noexcept { return m_r; }
  std::uint16_t layer() const noexcept { return m_layer; }
  float pixelClusterWidth() const noexcept { return m_pcw; }
  float localPositionY() const noexcept { return m_locPosY; }
  std::uint32_t sp_idx() const noexcept { return m_idx; }

  // writeable accessor
  float& x() noexcept { return m_x; }
  float& y() noexcept { return m_y; }
  float& phi() noexcept { return m_phi; }
  float& z() noexcept { return m_z; }
  float& r() noexcept { return m_r; }
  float& pixelClusterWidth() noexcept { return m_pcw; }
  float& localPositionY() noexcept { return m_locPosY; }
  std::uint32_t& sp_idx() noexcept { return m_idx; }

 private:
  float m_x{};
  float m_y{};
  float m_z{};
  float m_r{};
  float m_phi{};
  std::uint16_t m_layer{10000};
  std::uint32_t m_idx{std::numeric_limits<unsigned int>::max()};
  float m_pcw{};
  float m_locPosY{};
};

class GbtsEtaBin {
 public:
  struct CompareNodesByPhi {
    bool operator()(const GbtsNode* n1, const GbtsNode* n2) {
      return n1->phi() < n2->phi();
    }
  };

  GbtsEtaBin();

  void sortByPhi();
  void initializeNodes();
  bool empty() const { return m_vn.empty(); }

  void generatePhiIndexing(float dphi);

  float getMinBinRadius() const { return m_minRadius; }

  float getMaxBinRadius() const { return m_maxRadius; }

  /// nodes of the graph
  std::vector<const GbtsNode*> m_vn;
  std::vector<std::pair<float, unsigned int>> m_vPhiNodes;
  /// vectors of incoming edges, stores indices of edges in the edge vector
  std::vector<std::vector<unsigned int>> m_in;
  /// node attributes: m_minCutOnTau, m_maxCutOnTau, m_phi, m_r, m_z;
  std::vector<std::array<float, 5>> m_params;
  float m_minRadius{};
  float m_maxRadius{};

  unsigned int m_layerKey{0};
};

class GbtsDataStorage {
 public:
  explicit GbtsDataStorage(std::shared_ptr<const GbtsGeometry> geometry,
                           const SeedFinderGbtsConfig& config,
                           GbtsMLLookupTable mlLUT);

  int loadPixelGraphNodes(short layerIndex, std::span<const GbtsNode> coll,
                          bool useML);
  int loadStripGraphNodes(short layerIndex, std::span<const GbtsNode> coll);

  unsigned int numberOfNodes() const;
  void sortByPhi();
  void initializeNodes(bool useML);
  void generatePhiIndexing(float dphi);

  GbtsEtaBin& getEtaBin(int idx) {
    if (idx >= static_cast<int>(m_etaBins.size())) {
      idx = idx - 1;
    }
    return m_etaBins.at(idx);
  }

 protected:
  std::shared_ptr<const GbtsGeometry> m_geo;

  SeedFinderGbtsConfig m_config;

  GbtsMLLookupTable m_mlLUT;

  std::vector<GbtsEtaBin> m_etaBins;
};

class GbtsEdge {
 public:
  struct CompareLevel {
   public:
    bool operator()(const GbtsEdge* pE1, const GbtsEdge* pE2) {
      return pE1->m_level > pE2->m_level;
    }
  };

  GbtsEdge() = default;

  GbtsEdge(const GbtsNode* n1, const GbtsNode* n2, float p1, float p2, float p3)
      : m_n1(n1), m_n2(n2), m_level(1), m_next(1) {
    m_p[0] = p1;
    m_p[1] = p2;
    m_p[2] = p3;
  }

  const GbtsNode* m_n1{nullptr};
  const GbtsNode* m_n2{nullptr};

  std::int8_t m_level{-1};
  std::int8_t m_next{-1};

  std::uint8_t m_nNei{0};
  std::array<float, 3> m_p{};

  // global indices of the connected edges
  std::array<std::uint32_t, N_SEG_CONNS> m_vNei{};
};

}  // namespace Acts::Experimental
