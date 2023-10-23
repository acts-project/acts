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
#include "Acts/Seeding/GNN_Geometry.hpp"

#include <algorithm>
#include <map>
#include <vector>

namespace Acts {

constexpr size_t MAX_SEG_PER_NODE = 1000;  // was 30
constexpr size_t N_SEG_CONNS = 6;          // was 6

// new sp struct
template <typename space_point_t>
struct FTF_SP {
  const space_point_t *SP;  // want inside to have pointer
  int FTF_ID;
  int combined_ID;
  FTF_SP(const space_point_t *sp, int id, int combined_id)
      : SP(sp), FTF_ID(id), combined_ID{combined_id} {
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
class TrigFTF_GNN_Node {
 public:
  struct CompareByPhi {
    bool operator()(const TrigFTF_GNN_Node<space_point_t> *n1,
                    const TrigFTF_GNN_Node<space_point_t> *n2) {
      return (n1->m_sp_FTF.phi() < n2->m_sp_FTF.phi());
    }
  };

  TrigFTF_GNN_Node(const FTF_SP<space_point_t> &FTF_sp, float minT = -100.0,
                   float maxT = 100.0)
      : m_sp_FTF(FTF_sp), m_minCutOnTau(minT), m_maxCutOnTau(maxT) {}

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

  const FTF_SP<space_point_t> &m_sp_FTF;

  std::vector<unsigned int> m_in;  // indices of the edges in the edge storage
  std::vector<unsigned int> m_out;
  float m_minCutOnTau, m_maxCutOnTau;
};

template <typename space_point_t>
class TrigFTF_GNN_EtaBin {
 public:
  TrigFTF_GNN_EtaBin() { m_vn.clear(); }

  ~TrigFTF_GNN_EtaBin() {
    for (typename std::vector<TrigFTF_GNN_Node<space_point_t> *>::iterator it =
             m_vn.begin();
         it != m_vn.end(); ++it) {
      delete (*it);
    }
  }

  void sortByPhi() {
    std::sort(m_vn.begin(), m_vn.end(),
              typename Acts::TrigFTF_GNN_Node<space_point_t>::CompareByPhi());
  }

  bool empty() const { return m_vn.empty(); }

  void generatePhiIndexing(float dphi) {
    for (unsigned int nIdx = 0; nIdx < m_vn.size(); nIdx++) {
      TrigFTF_GNN_Node<space_point_t> *pN = m_vn.at(nIdx);
      // float phi = pN->m_sp.phi();
      // float phi = (std::atan(pN->m_sp.x() / pN->m_sp.y()));
      float phi = pN->m_sp_FTF.phi();
      if (phi <= M_PI - dphi) {
        continue;
      }

      m_vPhiNodes.push_back(
          std::pair<float, unsigned int>(phi - 2 * M_PI, nIdx));
    }

    for (unsigned int nIdx = 0; nIdx < m_vn.size(); nIdx++) {
      TrigFTF_GNN_Node<space_point_t> *pN = m_vn.at(nIdx);
      float phi = pN->m_sp_FTF.phi();
      m_vPhiNodes.push_back(std::pair<float, unsigned int>(phi, nIdx));
    }

    for (unsigned int nIdx = 0; nIdx < m_vn.size(); nIdx++) {
      TrigFTF_GNN_Node<space_point_t> *pN = m_vn.at(nIdx);
      float phi = pN->m_sp_FTF.phi();
      if (phi >= -M_PI + dphi) {
        break;
      }
      m_vPhiNodes.push_back(
          std::pair<float, unsigned int>(phi + 2 * M_PI, nIdx));
    }
  }

  std::vector<TrigFTF_GNN_Node<space_point_t> *> m_vn;
  // TODO change to
  // std::vector<std::unique_ptr<TrigFTF_GNN_Node<space_point_t>>> m_vn;
  std::vector<std::pair<float, unsigned int>> m_vPhiNodes;
};

template <typename space_point_t>
class TrigFTF_GNN_DataStorage {
 public:
  TrigFTF_GNN_DataStorage(const TrigFTF_GNN_Geometry<space_point_t> &g)
      : m_geo(g) {
    m_etaBins.reserve(g.num_bins());
    for (int k = 0; k < g.num_bins(); k++) {
      m_etaBins.emplace_back(TrigFTF_GNN_EtaBin<space_point_t>());
    }
  }

  int addSpacePoint(const FTF_SP<space_point_t> &sp, bool useClusterWidth) {
    const TrigFTF_GNN_Layer<space_point_t> *pL =
        m_geo.getTrigFTF_GNN_LayerByKey(sp.combined_ID);

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
        //   const Trk::SpacePoint* osp = sp.offlineSpacePoint();
        //   const InDet::PixelCluster* pCL = dynamic_cast<const
        //   InDet::PixelCluster*>(osp->clusterList().first);
        //   float cluster_width = pCL->width().widthPhiRZ().y();
        float cluster_width = 1;  // temporary while cluster width not available
        min_tau = 6.7 * (cluster_width - 0.2);
        max_tau =
            1.6 + 0.15 / (cluster_width + 0.2) + 6.1 * (cluster_width - 0.2);
      }

      m_etaBins.at(binIndex).m_vn.push_back(new TrigFTF_GNN_Node<space_point_t>(
          sp, min_tau, max_tau));  // adding ftf member to nodes
    } else {
      if (useClusterWidth) {
        //   const Trk::SpacePoint* osp = sp.offlineSpacePoint();
        //   const InDet::PixelCluster* pCL = dynamic_cast<const
        //   InDet::PixelCluster*>(osp->clusterList().first);
        //   float cluster_width = pCL->width().widthPhiRZ().y();
        float cluster_width = 1;  // temporary while cluster width not available
        if (cluster_width > 0.2) {
          return -3;
        }
      }
      m_etaBins.at(binIndex).m_vn.push_back(
          new TrigFTF_GNN_Node<space_point_t>(sp));
    }

    return 0;
  }

  // for safety to prevent passing as copy
  TrigFTF_GNN_DataStorage(const TrigFTF_GNN_DataStorage &) = delete;
  TrigFTF_GNN_DataStorage &operator=(const TrigFTF_GNN_DataStorage &) = delete;

  unsigned int numberOfNodes() const {
    unsigned int n = 0;

    for (auto &b : m_etaBins) {
      n += b.m_vn.size();
    }
    return n;
  }

  void getConnectingNodes(
      std::vector<const TrigFTF_GNN_Node<space_point_t> *> &vn) {
    vn.clear();
    vn.reserve(numberOfNodes());
    for (const auto &b : m_etaBins) {
      for (typename std::vector<
               TrigFTF_GNN_Node<space_point_t> *>::const_iterator nIt =
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

  const TrigFTF_GNN_EtaBin<space_point_t> &getEtaBin(int idx) const {
    if (idx >= static_cast<int>(m_etaBins.size())) {
      idx = idx - 1;
    }
    return m_etaBins.at(idx);
  }

 protected:
  const TrigFTF_GNN_Geometry<space_point_t> &m_geo;

  std::vector<TrigFTF_GNN_EtaBin<space_point_t>> m_etaBins;
};

template <typename space_point_t>
class TrigFTF_GNN_Edge {
 public:
  struct CompareLevel {
   public:
    bool operator()(const TrigFTF_GNN_Edge *pS1, const TrigFTF_GNN_Edge *pS2) {
      return pS1->m_level > pS2->m_level;
    }
  };

  TrigFTF_GNN_Edge(TrigFTF_GNN_Node<space_point_t> *n1,
                   TrigFTF_GNN_Node<space_point_t> *n2, float p1, float p2,
                   float p3, float p4)
      : m_n1(n1), m_n2(n2), m_level(1), m_next(1) {
    m_p[0] = p1;
    m_p[1] = p2;
    m_p[2] = p3;
    m_p[3] = p4;
  }

  TrigFTF_GNN_Edge() : m_n1(nullptr), m_n2(nullptr), m_level(-1), m_next(-1) {}

  // TrigFTF_GNN_Edge(const TrigFTF_GNN_Edge<space_point_t> &e)
  //     : m_n1(e.m_n1), m_n2(e.m_n2) {}

  // inline void initialize(TrigFTF_GNN_Node<space_point_t> *n1,
  //                        TrigFTF_GNN_Node<space_point_t> *n2) {
  //   m_n1 = n1;
  //   m_n2 = n2;
  //   m_level = 1;
  //   m_next = 1;
  //   m_nNei = 0;
  // }

  TrigFTF_GNN_Node<space_point_t> *m_n1{nullptr};
  TrigFTF_GNN_Node<space_point_t> *m_n2{nullptr};

  signed char m_level{}, m_next{};

  unsigned char m_nNei{0};
  float m_p[4]{};

  unsigned int m_vNei[N_SEG_CONNS]{};  // global indices of the connected edges
};

}  // namespace Acts
