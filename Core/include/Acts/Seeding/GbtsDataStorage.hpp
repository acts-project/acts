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
#include <limits>
#include <vector>
namespace Acts::Experimental {

#define MAX_SEG_PER_NODE 1000  // was 30
#define N_SEG_CONNS 6          // was 6

class GbtsGeometry;

struct GbtsNode {
  explicit GbtsNode(unsigned short l) : m_layer(l) {};

  inline float x() const { return m_x; }
  inline float y() const { return m_y; }

  inline float phi() const { return m_phi; }
  inline float z() const { return m_z; }
  inline float r() const { return m_r; }
  inline unsigned short layer() const { return m_layer; }
  inline float pixelClusterWidth() const { return m_pcw; }
  inline float localPositionY() const { return m_locPosY; }
  inline int sp_idx() const { return m_idx; }

  float m_x{}, m_y{}, m_z{}, m_r{}, m_phi{};
  unsigned short m_layer{10000};
  unsigned int m_idx{std::numeric_limits<unsigned int>::max()};
  float m_pcw{}, m_locPosY{};
};

class GbtsEtaBin {
 public:
  struct CompareNodesByPhi {
    bool operator()(const GbtsNode* n1, const GbtsNode* n2) {
      return n1->phi() < n2->phi();
    }
  };

  GbtsEtaBin();
  ~GbtsEtaBin();

  void sortByPhi();
  void initializeNodes();
  bool empty() const { return m_vn.empty(); }

  void generatePhiIndexing(float dphi);

  float getMinBinRadius() const { return m_minRadius; }

  float getMaxBinRadius() const { return m_maxRadius; }

  std::vector<const GbtsNode*> m_vn;  // nodes of the graph
  std::vector<std::pair<float, unsigned int>> m_vPhiNodes;
  std::vector<std::vector<unsigned int>>
      m_in;  // vectors of incoming edges, stores indices of edges in the edge
             // vector
  std::vector<std::array<float, 5>>
      m_params;  // node attributes: m_minCutOnTau, m_maxCutOnTau, m_phi, m_r,
                 // m_z;

  float m_minRadius{}, m_maxRadius{};

  unsigned int m_layerKey{0};
};

class GbtsDataStorage {
 public:
  explicit GbtsDataStorage(
      const GbtsGeometry& geometry, const SeedFinderGbtsConfig& config,
      const std::vector<std::array<float, 5>>& parsedLutFile);
  ~GbtsDataStorage();

  int loadPixelGraphNodes(short layerIndex, const std::vector<GbtsNode>& coll,
                          bool useML);
  int loadStripGraphNodes(short layerIndex, const std::vector<GbtsNode>& coll);

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
  const GbtsGeometry& m_geo;
  const SeedFinderGbtsConfig& m_config;
  const std::vector<std::array<float, 5>>& m_mlLUT;
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

  GbtsEdge(const GbtsNode* n1, const GbtsNode* n2, float p1, float p2, float p3)
      : m_n1(n1), m_n2(n2), m_level(1), m_next(1) {
    m_p[0] = p1;
    m_p[1] = p2;
    m_p[2] = p3;
  }

  GbtsEdge() = default;

  const GbtsNode* m_n1{nullptr};
  const GbtsNode* m_n2{nullptr};

  signed char m_level{-1}, m_next{-1};

  unsigned char m_nNei{0};
  float m_p[3]{};

  unsigned int m_vNei[N_SEG_CONNS]{};  // global indices of the connected edges
};

}  // namespace Acts::Experimental
