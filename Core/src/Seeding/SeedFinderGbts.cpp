// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SeedFinderGbts.hpp"

#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <numeric>
#include <utility>
#include <vector>

namespace Acts::Experimental {

SeedFinderGbts::SeedFinderGbts(
    SeedFinderGbtsConfig config, const GbtsGeometry* gbtsGeo,
    const std::vector<TrigInDetSiLayer>* layerGeometry,
    const GbtsLutParser* gbtsLutParser,
    std::unique_ptr<const Acts::Logger> logger)
    : m_config(std::move(config)),
      m_geo(gbtsGeo),
      m_layerGeometry(layerGeometry),
      m_lutParser(gbtsLutParser),
      m_logger(std::move(logger)) {
  m_config.phiSliceWidth = 2 * std::numbers::pi / m_config.nMaxPhiSlice;
}

SeedContainer2 SeedFinderGbts::CreateSeeds(
    const RoiDescriptor& roi,
    const SPContainerComponentsType& SpContainerComponents,
    int max_layers) const {
  auto& parsedLutFile = m_lutParser->getParsedLut();
  std::unique_ptr<GbtsDataStorage> storage =
      std::make_unique<GbtsDataStorage>(*m_geo, m_config, parsedLutFile);

  SeedContainer2 SeedContainer;
  std::vector<std::vector<GbtsNode>> node_storage =
      CreateNodes(SpContainerComponents, max_layers);
  unsigned int nPixelLoaded = 0;
  unsigned int nStripLoaded = 0;

  for (std::size_t l = 0; l < node_storage.size(); l++) {
    const std::vector<GbtsNode>& nodes = node_storage[l];

    if (nodes.empty()) {
      continue;
    }

    bool is_pixel = true;
    if (is_pixel) {  // placeholder for now until strip hits are added in

      nPixelLoaded += storage->loadPixelGraphNodes(l, nodes, m_config.useML);

    } else {
      nStripLoaded += storage->loadStripGraphNodes(l, nodes);
    }
  }
  ACTS_DEBUG("Loaded " << nPixelLoaded << " pixel spacepoints and "
                       << nStripLoaded << " strip spacepoints");

  storage->sortByPhi();

  storage->initializeNodes(m_config.useML);

  storage->generatePhiIndexing(1.5f * m_config.phiSliceWidth);

  std::vector<GbtsEdge> edgeStorage;

  std::pair<int, int> graphStats = buildTheGraph(roi, storage, edgeStorage);

  ACTS_DEBUG("Created graph with " << graphStats.first << " edges and "
                                   << graphStats.second << " edge links");

  if (graphStats.first == 0 || graphStats.second == 0) {
    ACTS_WARNING("Missing edges or edge connections");
  }

  int maxLevel = runCCA(graphStats.first, edgeStorage);

  ACTS_DEBUG("Reached Level " << maxLevel << " after GNN iterations");

  // std::vector<std::tuple<float, int, std::vector<unsigned int>>>
  // vSeedCandidates;
  std::vector<seedProperties> vSeedCandidates;
  extractSeedsFromTheGraph(maxLevel, graphStats.first,
                           std::get<0>(SpContainerComponents).size(),
                           edgeStorage, vSeedCandidates);

  if (vSeedCandidates.empty()) {
    ACTS_WARNING("No Seed Candidates");
  }

  for (const auto& seed : vSeedCandidates) {
    if (seed.isClone != 0) {
      continue;
    }

    // add to seed container:
    std::vector<SpacePointIndex2> Sp_Indexes{};
    Sp_Indexes.reserve(seed.spacepoints.size());

    for (const auto& sp_idx : seed.spacepoints) {
      Sp_Indexes.emplace_back(sp_idx);
    }

    auto newSeed = SeedContainer.createSeed();
    newSeed.assignSpacePointIndices(Sp_Indexes);
    newSeed.quality() = seed.seedQuality;
  }

  ACTS_DEBUG("GBTS created " << SeedContainer.size() << " seeds");

  return SeedContainer;
}

std::vector<std::vector<GbtsNode>> SeedFinderGbts::CreateNodes(
    const auto& container, int MaxLayers) const {
  std::vector<std::vector<GbtsNode>> node_storage(MaxLayers);
  // reserve for better efficiency

  for (auto& v : node_storage) {
    v.reserve(10000);
  }

  for (auto sp : std::get<0>(container)) {
    // for every sp in container,
    // add its variables to node_storage organised by layer
    int layer = sp.extra(std::get<1>(container));

    // add node to storage
    GbtsNode& node = node_storage[layer].emplace_back(layer);

    // fill the node with spacepoint variables

    node.m_x = sp.x();
    node.m_y = sp.y();
    node.m_z = sp.z();
    node.m_r = sp.r();
    node.m_phi = sp.phi();
    node.m_idx =
        sp.index();  // change node so that is uses SpacePointIndex2 (doesn't
                     // affect code but i guess it looks cleaner)
    node.m_pcw = sp.extra(std::get<2>(container));
    node.m_locPosY = sp.extra(std::get<3>(container));
  }

  return node_storage;
}

std::pair<int, int> SeedFinderGbts::buildTheGraph(
    const RoiDescriptor& roi, const std::unique_ptr<GbtsDataStorage>& storage,
    std::vector<GbtsEdge>& edgeStorage) const {
  // phi cut for triplets
  const float cut_dphi_max = m_config.LRTmode ? 0.07f : 0.012f;
  // curv cut for triplets
  const float cut_dcurv_max = m_config.LRTmode ? 0.015f : 0.001f;
  // tau cut for doublets and triplets
  const float cut_tau_ratio_max =
      m_config.LRTmode ? 0.015f : static_cast<float>(m_config.tau_ratio_cut);
  const float min_z0 = m_config.LRTmode ? -600.0 : roi.zedMinus();
  const float max_z0 = m_config.LRTmode ? 600.0 : roi.zedPlus();
  const float min_deltaPhi = m_config.LRTmode ? 0.01f : 0.001f;

  // used to calculate Z cut on doublets
  const float maxOuterRadius = m_config.LRTmode ? 1050.0 : 550.0;

  const float cut_zMinU = min_z0 + maxOuterRadius * roi.dzdrMinus();
  const float cut_zMaxU = max_z0 + maxOuterRadius * roi.dzdrPlus();

  // correction due to limited pT resolution
  float tripletPtMin = 0.8f * m_config.minPt;
  // to re-scale original tunings done for the 900 MeV pT cut
  const float pt_scale = 900.0f / m_config.minPt;

  float maxCurv = m_config.ptCoeff / tripletPtMin;

  float maxKappa_high_eta =
      m_config.LRTmode ? 1.0f * maxCurv : std::sqrt(0.8f) * maxCurv;
  float maxKappa_low_eta =
      m_config.LRTmode ? 1.0f * maxCurv : std::sqrt(0.6f) * maxCurv;

  // new settings for curvature cuts
  if (!m_config.useOldTunings && !m_config.LRTmode) {
    maxKappa_high_eta = 4.75e-4f * pt_scale;
    maxKappa_low_eta = 3.75e-4f * pt_scale;
  }

  const float dphi_coeff = m_config.LRTmode ? 1.0f * maxCurv : 0.68f * maxCurv;

  // the default sliding window along phi
  float deltaPhi = 0.5f * m_config.phiSliceWidth;

  unsigned int nConnections = 0;

  edgeStorage.reserve(m_config.nMaxEdges);

  int nEdges = 0;

  for (const auto& bg : m_geo->bin_groups()) {  // loop over bin groups

    GbtsEtaBin& B1 = storage->getEtaBin(bg.first);

    if (B1.empty()) {
      continue;
    }

    float rb1 = B1.getMinBinRadius();

    const unsigned int lk1 = B1.m_layerKey;

    for (const auto& b2_idx : bg.second) {
      const GbtsEtaBin& B2 = storage->getEtaBin(b2_idx);

      if (B2.empty()) {
        continue;
      }

      float rb2 = B2.getMaxBinRadius();

      if (m_config.useEtaBinning) {
        float abs_dr = std::fabs(rb2 - rb1);
        if (m_config.useOldTunings) {
          deltaPhi = min_deltaPhi + dphi_coeff * abs_dr;
        } else {
          if (abs_dr < 60.0) {
            deltaPhi = 0.002f + 4.33e-4f * pt_scale * abs_dr;
          } else {
            deltaPhi = 0.015f + 2.2e-4f * pt_scale * abs_dr;
          }
        }
      }

      unsigned int first_it = 0;

      for (unsigned int n1Idx = 0; n1Idx < B1.m_vn.size();
           n1Idx++) {  // loop over nodes in Layer 1

        std::vector<unsigned int>& v1In = B1.m_in[n1Idx];

        if (v1In.size() >= MAX_SEG_PER_NODE) {
          continue;
        }

        const std::array<float, 5>& n1pars = B1.m_params[n1Idx];

        float phi1 = n1pars[2];
        float r1 = n1pars[3];
        float z1 = n1pars[4];

        // sliding window phi1 +/- deltaPhi

        float minPhi = phi1 - deltaPhi;
        float maxPhi = phi1 + deltaPhi;

        for (unsigned int n2PhiIdx = first_it; n2PhiIdx < B2.m_vPhiNodes.size();
             n2PhiIdx++) {  // sliding window over nodes in Layer 2

          float phi2 = B2.m_vPhiNodes[n2PhiIdx].first;

          if (phi2 < minPhi) {
            first_it = n2PhiIdx;
            continue;
          }
          if (phi2 > maxPhi) {
            break;
          }

          unsigned int n2Idx = B2.m_vPhiNodes[n2PhiIdx].second;

          const std::vector<unsigned int>& v2In = B2.m_in[n2Idx];

          if (v2In.size() >= MAX_SEG_PER_NODE) {
            continue;
          }

          const std::array<float, 5>& n2pars = B2.m_params[n2Idx];

          float r2 = n2pars[3];

          float dr = r2 - r1;

          if (dr < m_config.minDeltaRadius) {
            continue;
          }

          float z2 = n2pars[4];

          float dz = z2 - z1;
          float tau = dz / dr;
          float ftau = std::fabs(tau);
          if (ftau > 36.0) {
            continue;
          }

          if (ftau < n1pars[0]) {
            continue;
          }
          if (ftau > n1pars[1]) {
            continue;
          }

          if (ftau < n2pars[0]) {
            continue;
          }
          if (ftau > n2pars[1]) {
            continue;
          }

          if (m_config.doubletFilterRZ) {
            float z0 = z1 - r1 * tau;

            if (z0 < min_z0 || z0 > max_z0) {
              continue;
            }

            float zouter = z0 + maxOuterRadius * tau;

            if (zouter < cut_zMinU || zouter > cut_zMaxU) {
              continue;
            }
          }

          float curv = (phi2 - phi1) / dr;
          float abs_curv = std::abs(curv);

          if (ftau < 4.0) {  // eta = 2.1
            if (abs_curv > maxKappa_low_eta) {
              continue;
            }
          } else {
            if (abs_curv > maxKappa_high_eta) {
              continue;
            }
          }

          float exp_eta = std::sqrt(1.f + tau * tau) - tau;

          if (m_config.matchBeforeCreate &&
              (lk1 == 80000 || lk1 == 81000)) {  // match edge candidate against
                                                 // edges incoming to n2

            bool isGood = v2In.size() <=
                          2;  // we must have enough incoming edges to decide

            if (!isGood) {
              float uat_1 = 1.0f / exp_eta;

              for (const auto& n2_in_idx : v2In) {
                float tau2 = edgeStorage.at(n2_in_idx).m_p[0];
                float tau_ratio = tau2 * uat_1 - 1.0f;

                if (std::abs(tau_ratio) >
                    m_config.tau_ratio_precut) {  // bad match
                  continue;
                }
                isGood = true;  // good match found
                break;
              }
            }

            if (!isGood) {  // no match found, skip creating [n1 <- n2] edge
              continue;
            }
          }

          float dPhi2 = curv * r2;
          float dPhi1 = curv * r1;

          if (nEdges < m_config.nMaxEdges) {
            edgeStorage.emplace_back(B1.m_vn[n1Idx], B2.m_vn[n2Idx], exp_eta,
                                     curv, phi1 + dPhi1);

            if (v1In.size() < MAX_SEG_PER_NODE) {
              v1In.push_back(nEdges);
            }

            int outEdgeIdx = nEdges;

            float uat_2 = 1.f / exp_eta;
            float Phi2 = phi2 + dPhi2;
            float curv2 = curv;

            for (const auto& inEdgeIdx :
                 v2In) {  // looking for neighbours of the new edge

              GbtsEdge* pS = &(edgeStorage.at(inEdgeIdx));

              if (pS->m_nNei >= N_SEG_CONNS) {
                continue;
              }

              float tau_ratio = pS->m_p[0] * uat_2 - 1.0f;

              if (std::abs(tau_ratio) > cut_tau_ratio_max) {  // bad match
                continue;
              }

              float dPhi = Phi2 - pS->m_p[2];

              if (dPhi < -std::numbers::pi) {
                dPhi += 2 * std::numbers::pi;
              } else if (dPhi > std::numbers::pi) {
                dPhi -= 2 * std::numbers::pi;
              }

              if (std::abs(dPhi) > cut_dphi_max) {
                continue;
              }

              float dcurv = curv2 - pS->m_p[1];

              if (dcurv < -cut_dcurv_max || dcurv > cut_dcurv_max) {
                continue;
              }

              pS->m_vNei[pS->m_nNei++] = outEdgeIdx;

              nConnections++;
            }
            nEdges++;
          }
        }  // loop over n2 (outer) nodes
      }  // loop over n1 (inner) nodes
    }  // loop over bins in Layer 2
  }  // loop over bin groups

  if (nEdges >= m_config.nMaxEdges) {
    ACTS_WARNING(
        "Maximum number of graph edges exceeded - possible efficiency loss "
        << nEdges);
  }
  return std::make_pair(nEdges, nConnections);
}

int SeedFinderGbts::runCCA(int nEdges,
                           std::vector<GbtsEdge>& edgeStorage) const {
  constexpr int maxIter = 15;

  int maxLevel = 0;

  int iter = 0;

  std::vector<GbtsEdge*> v_old;

  for (int edgeIndex = 0; edgeIndex < nEdges; edgeIndex++) {
    GbtsEdge* pS = &(edgeStorage[edgeIndex]);
    if (pS->m_nNei == 0) {
      continue;
    }

    v_old.push_back(pS);  // TO-DO: increment level for segments as they already
                          // have at least one neighbour
  }

  std::vector<GbtsEdge*> v_new;
  v_new.reserve(v_old.size());

  for (; iter < maxIter; iter++) {
    // generate proposals

    v_new.clear();

    for (auto pS : v_old) {
      int next_level = pS->m_level;

      for (int nIdx = 0; nIdx < pS->m_nNei; nIdx++) {
        unsigned int nextEdgeIdx = pS->m_vNei[nIdx];

        GbtsEdge* pN = &(edgeStorage[nextEdgeIdx]);

        if (pS->m_level == pN->m_level) {
          next_level = pS->m_level + 1;
          v_new.push_back(pS);
          break;
        }
      }

      pS->m_next = next_level;  // proposal
    }

    // update

    int nChanges = 0;

    for (auto pS : v_new) {
      if (pS->m_next != pS->m_level) {
        nChanges++;
        pS->m_level = pS->m_next;
        if (maxLevel < pS->m_level) {
          maxLevel = pS->m_level;
        }
      }
    }

    if (nChanges == 0) {
      break;
    }

    v_old.swap(v_new);
    v_new.clear();
  }

  return maxLevel;
}

void SeedFinderGbts::extractSeedsFromTheGraph(
    int maxLevel, int nEdges, int nHits, std::vector<GbtsEdge>& edgeStorage,
    std::vector<seedProperties>& vSeedCandidates) const {
  vSeedCandidates.clear();

  int minLevel = 3;  // a triplet + 2 confirmation

  if (m_config.LRTmode) {
    minLevel = 2;  // a triplet + 1 confirmation
  }

  if (maxLevel < minLevel) {
    return;
  }

  std::vector<GbtsEdge*> vSeeds;

  vSeeds.reserve(nEdges / 2);

  for (int edgeIndex = 0; edgeIndex < nEdges; edgeIndex++) {
    GbtsEdge* pS = &(edgeStorage.at(edgeIndex));

    if (pS->m_level < minLevel) {
      continue;
    }

    vSeeds.push_back(pS);
  }

  if (vSeeds.empty()) {
    return;
  }

  std::sort(vSeeds.begin(), vSeeds.end(), GbtsEdge::CompareLevel());

  // backtracking

  vSeedCandidates.reserve(vSeeds.size());

  GbtsTrackingFilter tFilter(*m_layerGeometry, edgeStorage, m_config);

  for (auto pS : vSeeds) {
    if (pS->m_level == -1) {
      continue;
    }

    GbtsEdgeState rs(false);

    tFilter.followTrack(pS, rs);

    if (!rs.m_initialized) {
      continue;
    }

    if (static_cast<int>(rs.m_vs.size()) < minLevel) {
      continue;
    }

    float seed_eta = std::abs(-std::log(pS->m_p[0]));

    std::vector<const GbtsNode*> vN;

    for (std::vector<GbtsEdge*>::reverse_iterator sIt = rs.m_vs.rbegin();
         sIt != rs.m_vs.rend(); ++sIt) {
      if (seed_eta > m_config.edge_mask_min_eta) {
        (*sIt)->m_level = -1;  // mark as collected
      }

      if (sIt == rs.m_vs.rbegin()) {
        vN.push_back((*sIt)->m_n1);
      }

      vN.push_back((*sIt)->m_n2);
    }

    if (vN.size() < 3) {
      continue;
    }

    std::vector<unsigned int> vSpIdx;

    vSpIdx.resize(vN.size());

    for (unsigned int k = 0; k < vN.size(); k++) {
      vSpIdx[k] = vN[k]->sp_idx();
    }

    vSeedCandidates.emplace_back(-rs.m_J / vN.size(), 0, std::move(vSpIdx));
  }

  // clone removal code goes below ...

  std::sort(vSeedCandidates.begin(), vSeedCandidates.end());

  std::vector<int> vTrackIds(vSeedCandidates.size());

  // fills the vector from 1 to N

  std::iota(vTrackIds.begin(), vTrackIds.end(), 1);

  std::vector<int> H2T(nHits + 1, 0);  // hit to track associations

  int seedIdx = 0;

  for (const auto& seed : vSeedCandidates) {
    for (const auto& h : seed.spacepoints) {  // loop over spacepoints indices

      unsigned int hit_id = h + 1;

      int tid = H2T[hit_id];
      int trackId = vTrackIds[seedIdx];

      if (tid == 0 || tid > trackId) {  // un-used hit or used by a lesser track

        H2T[hit_id] = trackId;  // overwrite
      }
    }

    seedIdx++;
  }

  for (unsigned int trackIdx = 0; trackIdx < vSeedCandidates.size();
       trackIdx++) {
    int nTotal = vSeedCandidates[trackIdx].spacepoints.size();
    int nOther = 0;

    int trackId = vTrackIds[trackIdx];

    for (const auto& h : vSeedCandidates[trackIdx].spacepoints) {
      unsigned int hit_id = h + 1;

      int tid = H2T[hit_id];

      if (tid != trackId) {  // taken by a better candidate
        nOther++;
      }
    }

    if (nOther > m_config.hit_share_threshold * nTotal) {
      vSeedCandidates[trackIdx].isClone = -1;  // reject
    }
  }
}

}  // namespace Acts::Experimental
