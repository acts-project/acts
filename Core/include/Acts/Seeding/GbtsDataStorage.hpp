// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"

#include <algorithm>
#include <map>
#include <vector>

namespace Acts {

constexpr std::size_t MAX_SEG_PER_NODE = 1000;  // was 30
constexpr std::size_t N_SEG_CONNS = 6;          // was 6

// new sp struct
template <typename space_point_t>
struct GbtsSP {
  const space_point_t *SP;  // want inside to have pointer
  int gbtsID;
  int combined_ID;
  GbtsSP(const space_point_t *sp, int id, int combined_id)
      : SP(sp), gbtsID(id), combined_ID{combined_id} {
    if (SP->sourceLinks().size() == 1) {  // pixels have 1 SL
      m_isPixel = true;
    } else {
      m_isPixel = false;
    }
    m_phi = std::atan(SP->x() / SP->y());
  };
  bool isPixel() const { return m_isPixel; }
  bool isSCT() const { return !m_isPixel; }
  float phi() const { return m_phi; }
  bool m_isPixel;
  float m_phi;
};

template <typename space_point_t>
class GbtsNode {
 public:
  struct CompareByPhi {
    bool operator()(const GbtsNode<space_point_t> *n1,
                    const GbtsNode<space_point_t> *n2) {
      return (n1->m_spGbts.phi() < n2->m_spGbts.phi());
    }
  };

  GbtsNode(const GbtsSP<space_point_t> &spGbts, float minT = -100.0,
           float maxT = 100.0)
      : m_spGbts(spGbts), m_minCutOnTau(minT), m_maxCutOnTau(maxT) {}

  inline void addIn(int i) {
    if (m_in.size() < MAX_SEG_PER_NODE) {
      m_in.push_back(i);
    }
  }

  inline void addOut(int i) {
    if (m_out.size() < MAX_SEG_PER_NODE) {
      m_out.push_back(i);
    }
  }

  inline bool isConnector() const {
    if (m_in.empty() || m_out.empty()) {
      return false;
    }
    return true;
  }

  inline bool isFull() const {
    if (m_in.size() == MAX_SEG_PER_NODE && m_out.size() == MAX_SEG_PER_NODE) {
      return true;
    } else {
      return false;
    }
  }

  const GbtsSP<space_point_t> &m_spGbts;

  std::vector<unsigned int> m_in;  // indices of the edges in the edge storage
  std::vector<unsigned int> m_out;
  float m_minCutOnTau, m_maxCutOnTau;
};

template <typename space_point_t>
class GbtsEtaBin {
 public:
  GbtsEtaBin() { m_vn.clear(); }

  ~GbtsEtaBin() {
    for (typename std::vector<GbtsNode<space_point_t> *>::iterator it =
             m_vn.begin();
         it != m_vn.end(); ++it) {
      delete (*it);
    }
  }

  void sortByPhi() {
    std::sort(m_vn.begin(), m_vn.end(),
              typename Acts::GbtsNode<space_point_t>::CompareByPhi());
  }

  bool empty() const { return m_vn.empty(); }

  void generatePhiIndexing(float dphi) {
    for (unsigned int nIdx = 0; nIdx < m_vn.size(); nIdx++) {
      GbtsNode<space_point_t> *pN = m_vn.at(nIdx);
      // float phi = pN->m_sp.phi();
      // float phi = (std::atan(pN->m_sp.x() / pN->m_sp.y()));
      float phi = pN->m_spGbts.phi();
      if (phi <= M_PI - dphi) {
        continue;
      }

      m_vPhiNodes.push_back(
          std::pair<float, unsigned int>(phi - 2 * M_PI, nIdx));
    }

    for (unsigned int nIdx = 0; nIdx < m_vn.size(); nIdx++) {
      GbtsNode<space_point_t> *pN = m_vn.at(nIdx);
      float phi = pN->m_spGbts.phi();
      m_vPhiNodes.push_back(std::pair<float, unsigned int>(phi, nIdx));
    }

    for (unsigned int nIdx = 0; nIdx < m_vn.size(); nIdx++) {
      GbtsNode<space_point_t> *pN = m_vn.at(nIdx);
      float phi = pN->m_spGbts.phi();
      if (phi >= -M_PI + dphi) {
        break;
      }
      m_vPhiNodes.push_back(
          std::pair<float, unsigned int>(phi + 2 * M_PI, nIdx));
    }
  }

  std::vector<GbtsNode<space_point_t> *> m_vn;
  // TODO change to
  // std::vector<std::unique_ptr<GbtsNode<space_point_t>>> m_vn;
  std::vector<std::pair<float, unsigned int>> m_vPhiNodes;
};

template <typename space_point_t>
class GbtsDataStorage {
 public:
  GbtsDataStorage(const GbtsGeometry<space_point_t> &g) : m_geo(g) {
    m_etaBins.reserve(g.num_bins());
    for (int k = 0; k < g.num_bins(); k++) {
      m_etaBins.emplace_back(GbtsEtaBin<space_point_t>());
    }
  }

  int addSpacePoint(const GbtsSP<space_point_t> &sp, bool useClusterWidth) {
    const GbtsLayer<space_point_t> *pL =
        m_geo.getGbtsLayerByKey(sp.combined_ID);

    if (pL == nullptr) {
      return -1;
    }

    int binIndex = pL->getEtaBin(sp.SP->z(), sp.SP->r());

    if (binIndex == -1) {
      return -2;
    }

    bool isBarrel = (pL->m_layer.m_type == 0);

    if (isBarrel) {
      float min_tau = -100.0;
      float max_tau = 100.0;
      // can't do this bit yet as dont have cluster width
      if (useClusterWidth) {
        float cluster_width = 1;  // temporary while cluster width not available
        min_tau = 6.7 * (cluster_width - 0.2);
        max_tau =
            1.6 + 0.15 / (cluster_width + 0.2) + 6.1 * (cluster_width - 0.2);
      }

      m_etaBins.at(binIndex).m_vn.push_back(new GbtsNode<space_point_t>(
          sp, min_tau, max_tau));  // adding ftf member to nodes
    } else {
      if (useClusterWidth) {
        float cluster_width = 1;  // temporary while cluster width not available
        if (cluster_width > 0.2) {
          return -3;
        }
      }
      m_etaBins.at(binIndex).m_vn.push_back(new GbtsNode<space_point_t>(sp));
    }

    return 0;
  }

  // for safety to prevent passing as copy
  GbtsDataStorage(const GbtsDataStorage &) = delete;
  GbtsDataStorage &operator=(const GbtsDataStorage &) = delete;

  unsigned int numberOfNodes() const {
    unsigned int n = 0;

    for (auto &b : m_etaBins) {
      n += b.m_vn.size();
    }
    return n;
  }

  void getConnectingNodes(std::vector<const GbtsNode<space_point_t> *> &vn) {
    vn.clear();
    vn.reserve(numberOfNodes());
    for (const auto &b : m_etaBins) {
      for (typename std::vector<GbtsNode<space_point_t> *>::const_iterator nIt =
               b.m_vn.begin();
           nIt != b.m_vn.end(); ++nIt) {
        if ((*nIt)->m_in.empty()) {
          continue;
        }
        if ((*nIt)->m_out.empty()) {
          continue;
        }
        vn.push_back(*nIt);
      }
    }
  }

  void sortByPhi() {
    for (auto &b : m_etaBins) {
      b.sortByPhi();
    }
  }

  void generatePhiIndexing(float dphi) {
    for (auto &b : m_etaBins) {
      b.generatePhiIndexing(dphi);
    }
  }

  const GbtsEtaBin<space_point_t> &getEtaBin(int idx) const {
    if (idx >= static_cast<int>(m_etaBins.size())) {
      idx = idx - 1;
    }
    return m_etaBins.at(idx);
  }

 protected:
  const GbtsGeometry<space_point_t> &m_geo;

  std::vector<GbtsEtaBin<space_point_t>> m_etaBins;
};

template <typename space_point_t>
class GbtsEdge {
 public:
  struct CompareLevel {
   public:
    bool operator()(const GbtsEdge *pS1, const GbtsEdge *pS2) {
      return pS1->m_level > pS2->m_level;
    }
  };

  GbtsEdge(GbtsNode<space_point_t> *n1, GbtsNode<space_point_t> *n2, float p1,
           float p2, float p3, float p4)
      : m_n1(n1), m_n2(n2), m_level(1), m_next(1) {
    m_p[0] = p1;
    m_p[1] = p2;
    m_p[2] = p3;
    m_p[3] = p4;
  }

  GbtsEdge() : m_n1(nullptr), m_n2(nullptr), m_level(-1), m_next(-1) {}

  // GbtsEdge(const GbtsEdge<space_point_t> &e)
  //     : m_n1(e.m_n1), m_n2(e.m_n2) {}

  // inline void initialize(GbtsNode<space_point_t> *n1,
  //                        GbtsNode<space_point_t> *n2) {
  //   m_n1 = n1;
  //   m_n2 = n2;
  //   m_level = 1;
  //   m_next = 1;
  //   m_nNei = 0;
  // }

  GbtsNode<space_point_t> *m_n1{nullptr};
  GbtsNode<space_point_t> *m_n2{nullptr};

  signed char m_level{}, m_next{};

  unsigned char m_nNei{0};
  float m_p[4]{};

  unsigned int m_vNei[N_SEG_CONNS]{};  // global indices of the connected edges
};

}  // namespace Acts
