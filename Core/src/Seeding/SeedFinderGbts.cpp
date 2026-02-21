// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SeedFinderGbts.hpp"

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <numbers>
#include <numeric>
#include <utility>
#include <vector>

namespace Acts::Experimental {

SeedFinderGbts::SeedFinderGbts(
    const GbtsConfig& config, std::unique_ptr<GbtsGeometry> gbtsGeo,
    const std::vector<TrigInDetSiLayer>& layerGeometry,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(config),
      m_geo(std::move(gbtsGeo)),
      m_layerGeometry(&layerGeometry),
      m_logger(std::move(logger)) {
  m_cfg.phiSliceWidth = 2 * std::numbers::pi_v<float> / m_cfg.nMaxPhiSlice;

  m_mlLut = parseGbtsMlLookupTable(m_cfg.lutInputFile);
}

SeedContainer2 SeedFinderGbts::createSeeds(
    const RoiDescriptor& roi, const SpacePointContainer2& spacePoints,
    const std::uint32_t maxLayers) const {
  auto storage = std::make_unique<GbtsDataStorage>(m_cfg, m_geo, m_mlLut);

  SeedContainer2 SeedContainer;
  std::vector<std::vector<GbtsNode>> node_storage =
      createNodes(spacePoints, maxLayers);
  std::uint32_t nPixelLoaded = 0;
  std::uint32_t nStripLoaded = 0;

  for (std::uint16_t l = 0; l < node_storage.size(); l++) {
    const std::vector<GbtsNode>& nodes = node_storage[l];

    if (nodes.empty()) {
      continue;
    }

    bool is_pixel = true;
    // placeholder for now until strip hits are added in
    if (is_pixel) {
      nPixelLoaded += storage->loadPixelGraphNodes(l, nodes, m_cfg.useMl);
    } else {
      nStripLoaded += storage->loadStripGraphNodes(l, nodes);
    }
  }
  ACTS_DEBUG("Loaded " << nPixelLoaded << " pixel space points and "
                       << nStripLoaded << " strip space points");

  storage->sortByPhi();

  storage->initializeNodes(m_cfg.useMl);

  storage->generatePhiIndexing(1.5f * m_cfg.phiSliceWidth);

  std::vector<GbtsEdge> edgeStorage;

  std::pair<std::int32_t, std::int32_t> graphStats =
      buildTheGraph(roi, storage, edgeStorage);

  ACTS_DEBUG("Created graph with " << graphStats.first << " edges and "
                                   << graphStats.second << " edge links");

  if (graphStats.first == 0 || graphStats.second == 0) {
    ACTS_WARNING("Missing edges or edge connections");
  }

  std::uint32_t maxLevel = runCCA(graphStats.first, edgeStorage);

  ACTS_DEBUG("Reached Level " << maxLevel << " after GNN iterations");

  std::vector<SeedProperties> vSeedCandidates;
  extractSeedsFromTheGraph(maxLevel, graphStats.first, spacePoints.size(),
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
    Sp_Indexes.reserve(seed.spacePoints.size());

    for (const auto& sp_idx : seed.spacePoints) {
      Sp_Indexes.emplace_back(sp_idx);
    }

    auto newSeed = SeedContainer.createSeed();
    newSeed.assignSpacePointIndices(Sp_Indexes);
    newSeed.quality() = seed.seedQuality;
  }

  ACTS_DEBUG("GBTS created " << SeedContainer.size() << " seeds");

  return SeedContainer;
}

GbtsMlLookupTable SeedFinderGbts::parseGbtsMlLookupTable(
    const std::string& lutInputFile) {
  if (!m_cfg.useMl) {
    return {};
  }
  if (lutInputFile.empty()) {
    throw std::runtime_error("Cannot find ML predictor LUT file");
  }

  std::ifstream ifs(std::string(lutInputFile).c_str());
  if (!ifs.is_open()) {
    throw std::runtime_error("Failed to open LUT file");
  }

  GbtsMlLookupTable mlLut;
  mlLut.reserve(100);

  float clWidth{};
  float min1{};
  float max1{};
  float min2{};
  float max2{};
  while (ifs >> clWidth >> min1 >> max1 >> min2 >> max2) {
    mlLut.emplace_back(std::array<float, 5>{clWidth, min1, max1, min2, max2});
  }

  if (!ifs.eof()) {
    // ended if parse error present, not clean EOF
    throw std::runtime_error("Stopped reading LUT file due to parse error");
  }

  return mlLut;
}

std::vector<std::vector<GbtsNode>> SeedFinderGbts::createNodes(
    const SpacePointContainer2& spacePoints,
    const std::uint32_t maxLayers) const {
  auto layerColumn = spacePoints.column<std::uint32_t>("layerId");
  auto clusterWidthColomn = spacePoints.column<float>("clusterWidth");
  auto localPositionColomn = spacePoints.column<float>("localPositionY");

  std::vector<std::vector<GbtsNode>> node_storage(maxLayers);
  // reserve for better efficiency
  for (auto& v : node_storage) {
    v.reserve(10000);
  }

  for (auto sp : spacePoints) {
    // for every sp in container,
    // add its variables to node_storage organised by layer
    std::uint16_t layer = sp.extra(layerColumn);

    // add node to storage
    GbtsNode& node = node_storage[layer].emplace_back(layer);

    // fill the node with space point variables

    node.x = sp.x();
    node.y = sp.y();
    node.z = sp.z();
    node.r = sp.r();
    node.phi = sp.phi();
    node.idx = sp.index();
    node.pcw = sp.extra(clusterWidthColomn);
    node.locPosY = sp.extra(localPositionColomn);
  }

  return node_storage;
}

std::pair<std::int32_t, std::int32_t> SeedFinderGbts::buildTheGraph(
    const RoiDescriptor& roi, const std::unique_ptr<GbtsDataStorage>& storage,
    std::vector<GbtsEdge>& edgeStorage) const {
  // phi cut for triplets
  const float cutDPhiMax = m_cfg.lrtMode ? 0.07f : 0.012f;
  // curv cut for triplets
  const float cutDCurvMax = m_cfg.lrtMode ? 0.015f : 0.001f;
  // tau cut for doublets and triplets
  const float cutTauRatioMax = m_cfg.lrtMode ? 0.015f : m_cfg.tauRatioCut;
  const float minZ0 =
      m_cfg.lrtMode ? -600.0f : static_cast<float>(roi.zedMinus());
  const float maxZ0 =
      m_cfg.lrtMode ? 600.0f : static_cast<float>(roi.zedPlus());
  const float minDeltaPhi = m_cfg.lrtMode ? 0.01f : 0.001f;

  // used to calculate Z cut on doublets
  const float maxOuterRadius = m_cfg.lrtMode ? 1050.0f : 550.0f;

  const float cutZMinU =
      minZ0 + maxOuterRadius * static_cast<float>(roi.dzdrMinus());
  const float cutZMaxU =
      maxZ0 + maxOuterRadius * static_cast<float>(roi.dzdrPlus());

  // correction due to limited pT resolution
  float tripletPtMin = 0.8f * m_cfg.minPt;

  // to re-scale original tunings done for the 900 MeV pT cut
  const float ptScale = (0.9f * UnitConstants::GeV) / m_cfg.minPt;

  const float maxCurv = m_cfg.ptCoeff / tripletPtMin;

  float maxKappaHighEta =
      m_cfg.lrtMode ? 1.0f * maxCurv : std::sqrt(0.8f) * maxCurv;
  float maxKappaLowEta =
      m_cfg.lrtMode ? 1.0f * maxCurv : std::sqrt(0.6f) * maxCurv;

  // new settings for curvature cuts
  if (!m_cfg.useOldTunings && !m_cfg.lrtMode) {
    maxKappaHighEta = 4.75e-4f * ptScale;
    maxKappaLowEta = 3.75e-4f * ptScale;
  }

  const float dPhiCoeff = m_cfg.lrtMode ? 1.0f * maxCurv : 0.68f * maxCurv;

  // the default sliding window along phi
  float deltaPhi = 0.5f * m_cfg.phiSliceWidth;

  std::uint32_t nConnections = 0;

  edgeStorage.reserve(m_cfg.nMaxEdges);

  std::uint32_t nEdges = 0;

  // loop over bin groups
  for (const auto& bg : m_geo->binGroups()) {
    GbtsEtaBin& B1 = storage->getEtaBin(bg.first);

    if (B1.empty()) {
      continue;
    }

    const float rb1 = B1.minRadius;

    const std::uint32_t lk1 = B1.layerKey;

    for (const int b2Idx : bg.second) {
      const GbtsEtaBin& B2 = storage->getEtaBin(b2Idx);

      if (B2.empty()) {
        continue;
      }

      const float rb2 = B2.maxRadius;

      if (m_cfg.useEtaBinning) {
        const float absDr = std::fabs(rb2 - rb1);
        if (m_cfg.useOldTunings) {
          deltaPhi = minDeltaPhi + dPhiCoeff * absDr;
        } else {
          if (absDr < 60.0) {
            deltaPhi = 0.002f + 4.33e-4f * ptScale * absDr;
          } else {
            deltaPhi = 0.015f + 2.2e-4f * ptScale * absDr;
          }
        }
      }

      std::uint32_t firstIt = 0;
      // loop over nodes in Layer 1
      for (std::uint32_t n1Idx = 0; n1Idx < B1.vn.size(); ++n1Idx) {
        std::vector<std::uint32_t>& v1In = B1.in[n1Idx];

        if (v1In.size() >= gbtsMaxSegPerNode) {
          continue;
        }

        const std::array<float, 5>& n1pars = B1.params[n1Idx];

        const float phi1 = n1pars[2];
        const float r1 = n1pars[3];
        const float z1 = n1pars[4];

        // sliding window phi1 +/- deltaPhi

        const float minPhi = phi1 - deltaPhi;
        const float maxPhi = phi1 + deltaPhi;

        // sliding window over nodes in Layer 2
        for (std::uint32_t n2PhiIdx = firstIt; n2PhiIdx < B2.vPhiNodes.size();
             n2PhiIdx++) {
          const float phi2 = B2.vPhiNodes[n2PhiIdx].first;

          if (phi2 < minPhi) {
            firstIt = n2PhiIdx;
            continue;
          }
          if (phi2 > maxPhi) {
            break;
          }

          const std::uint32_t n2Idx = B2.vPhiNodes[n2PhiIdx].second;

          const std::vector<std::uint32_t>& v2In = B2.in[n2Idx];

          if (v2In.size() >= gbtsMaxSegPerNode) {
            continue;
          }

          const std::array<float, 5>& n2pars = B2.params[n2Idx];

          const float r2 = n2pars[3];

          const float dr = r2 - r1;

          if (dr < m_cfg.minDeltaRadius) {
            continue;
          }

          const float z2 = n2pars[4];

          const float dz = z2 - z1;
          const float tau = dz / dr;
          const float ftau = std::fabs(tau);
          if (ftau > 36.0f) {
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

          if (m_cfg.doubletFilterRZ) {
            const float z0 = z1 - r1 * tau;

            if (z0 < minZ0 || z0 > maxZ0) {
              continue;
            }

            const float zouter = z0 + maxOuterRadius * tau;

            if (zouter < cutZMinU || zouter > cutZMaxU) {
              continue;
            }
          }

          const float curv = (phi2 - phi1) / dr;
          const float absCurv = std::abs(curv);

          if (ftau < 4.0f) {  // eta = 2.1
            if (absCurv > maxKappaLowEta) {
              continue;
            }
          } else {
            if (absCurv > maxKappaHighEta) {
              continue;
            }
          }

          const float expEta = std::sqrt(1.f + tau * tau) - tau;

          // match edge candidate against edges incoming to n2
          if (m_cfg.matchBeforeCreate && (lk1 == 80000 || lk1 == 81000)) {
            // we must have enough incoming edges to decide
            bool isGood = v2In.size() <= 2;

            if (!isGood) {
              const float uat1 = 1.0f / expEta;

              for (const auto& n2InIdx : v2In) {
                const float tau2 = edgeStorage.at(n2InIdx).p[0];
                const float tauRatio = tau2 * uat1 - 1.0f;

                if (std::abs(tauRatio) > m_cfg.tauRatioPrecut) {  // bad match
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

          const float dPhi2 = curv * r2;
          const float dPhi1 = curv * r1;

          if (nEdges < m_cfg.nMaxEdges) {
            edgeStorage.emplace_back(B1.vn[n1Idx], B2.vn[n2Idx], expEta, curv,
                                     phi1 + dPhi1);

            if (v1In.size() < gbtsMaxSegPerNode) {
              v1In.push_back(nEdges);
            }

            const std::uint32_t outEdgeIdx = nEdges;

            const float uat2 = 1.f / expEta;
            const float phi2u = phi2 + dPhi2;
            const float curv2 = curv;

            // looking for neighbours of the new edge
            for (const auto& inEdgeIdx : v2In) {
              GbtsEdge* pS = &(edgeStorage.at(inEdgeIdx));

              if (pS->nNei >= gbtsNumSegConns) {
                continue;
              }

              const float tauRatio = pS->p[0] * uat2 - 1.0f;

              if (std::abs(tauRatio) > cutTauRatioMax) {  // bad match
                continue;
              }

              float dPhi = phi2u - pS->p[2];

              if (dPhi < -std::numbers::pi_v<float>) {
                dPhi += 2 * std::numbers::pi_v<float>;
              } else if (dPhi > std::numbers::pi_v<float>) {
                dPhi -= 2 * std::numbers::pi_v<float>;
              }

              if (std::abs(dPhi) > cutDPhiMax) {
                continue;
              }

              const float dcurv = curv2 - pS->p[1];

              if (dcurv < -cutDCurvMax || dcurv > cutDCurvMax) {
                continue;
              }

              pS->vNei[pS->nNei] = outEdgeIdx;
              ++pS->nNei;

              nConnections++;
            }
            nEdges++;
          }
        }  // loop over n2 (outer) nodes
      }  // loop over n1 (inner) nodes
    }  // loop over bins in Layer 2
  }  // loop over bin groups

  if (nEdges >= m_cfg.nMaxEdges) {
    ACTS_WARNING(
        "Maximum number of graph edges exceeded - possible efficiency loss "
        << nEdges);
  }
  return std::make_pair(nEdges, nConnections);
}

std::int32_t SeedFinderGbts::runCCA(const std::uint32_t nEdges,
                                    std::vector<GbtsEdge>& edgeStorage) const {
  constexpr std::uint32_t maxIter = 15;

  std::int32_t maxLevel = 0;

  std::uint32_t iter = 0;

  std::vector<GbtsEdge*> v_old;

  for (std::uint32_t edgeIndex = 0; edgeIndex < nEdges; ++edgeIndex) {
    GbtsEdge* pS = &(edgeStorage[edgeIndex]);
    if (pS->nNei == 0) {
      continue;
    }

    // TODO: increment level for segments as they already have at least one
    // neighbour
    v_old.push_back(pS);
  }

  std::vector<GbtsEdge*> v_new;
  v_new.reserve(v_old.size());

  // generate proposals
  for (; iter < maxIter; iter++) {
    v_new.clear();

    for (GbtsEdge* pS : v_old) {
      std::int32_t nextLevel = pS->level;

      for (std::uint32_t nIdx = 0; nIdx < pS->nNei; ++nIdx) {
        std::uint32_t nextEdgeIdx = pS->vNei[nIdx];

        GbtsEdge* pN = &(edgeStorage[nextEdgeIdx]);

        if (pS->level == pN->level) {
          nextLevel = pS->level + 1;
          v_new.push_back(pS);
          break;
        }
      }

      // proposal
      pS->next = static_cast<std::int8_t>(nextLevel);
    }

    // update

    std::uint32_t nChanges = 0;

    for (auto pS : v_new) {
      if (pS->next != pS->level) {
        nChanges++;
        pS->level = pS->next;
        if (maxLevel < pS->level) {
          maxLevel = pS->level;
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
    std::uint32_t maxLevel, std::uint32_t nEdges, std::int32_t nHits,
    std::vector<GbtsEdge>& edgeStorage,
    std::vector<SeedProperties>& vSeedCandidates) const {
  vSeedCandidates.clear();

  // a triplet + 2 confirmation
  std::uint32_t minLevel = 3;

  if (m_cfg.lrtMode) {
    // a triplet + 1 confirmation
    minLevel = 2;
  }

  if (maxLevel < minLevel) {
    return;
  }

  std::vector<GbtsEdge*> vSeeds;

  vSeeds.reserve(nEdges / 2);

  for (std::uint32_t edgeIndex = 0; edgeIndex < nEdges; ++edgeIndex) {
    GbtsEdge* pS = &(edgeStorage.at(edgeIndex));

    if (pS->level < static_cast<std::int32_t>(minLevel)) {
      continue;
    }

    vSeeds.push_back(pS);
  }

  if (vSeeds.empty()) {
    return;
  }

  std::ranges::sort(vSeeds, std::ranges::greater{},
                    [](const GbtsEdge* e) { return e->level; });

  // backtracking

  vSeedCandidates.reserve(vSeeds.size());

  GbtsTrackingFilter tFilter(m_cfg, *m_layerGeometry, edgeStorage);

  for (GbtsEdge* pS : vSeeds) {
    if (pS->level == -1) {
      continue;
    }

    GbtsEdgeState rs = tFilter.followTrack(*pS);

    if (!rs.initialized) {
      continue;
    }

    if (minLevel > static_cast<std::uint32_t>(rs.vs.size())) {
      continue;
    }

    const float seedEta = std::abs(-std::log(pS->p[0]));

    std::vector<const GbtsNode*> vN;

    for (auto sIt = rs.vs.rbegin(); sIt != rs.vs.rend(); ++sIt) {
      if (seedEta > m_cfg.edgeMaskMinEta) {
        // mark as collected
        (*sIt)->level = -1;
      }

      if (sIt == rs.vs.rbegin()) {
        vN.push_back((*sIt)->n1);
      }

      vN.push_back((*sIt)->n2);
    }

    if (vN.size() < 3) {
      continue;
    }

    std::vector<std::uint32_t> vSpIdx;
    vSpIdx.resize(vN.size());

    for (std::uint32_t k = 0; k < vN.size(); ++k) {
      vSpIdx[k] = vN[k]->idx;
    }

    vSeedCandidates.emplace_back(-rs.j / vN.size(), 0, std::move(vSpIdx));
  }

  // clone removal code goes below ...

  std::ranges::sort(vSeedCandidates,
                    [](const SeedProperties& s1, const SeedProperties& s2) {
                      return s1 < s2;
                    });

  std::vector<std::int32_t> vTrackIds(vSeedCandidates.size());

  // fills the vector from 1 to N

  std::iota(vTrackIds.begin(), vTrackIds.end(), 1);

  std::vector<std::uint32_t> h2t(nHits + 1, 0);  // hit to track associations

  std::uint32_t seedIdx = 0;

  for (const auto& seed : vSeedCandidates) {
    // loop over space points indices
    for (const auto& h : seed.spacePoints) {
      const std::uint32_t hitId = h + 1;

      const std::uint32_t tid = h2t[hitId];
      const std::uint32_t trackId = vTrackIds[seedIdx];

      // un-used hit or used by a lesser track
      if (tid == 0 || tid > trackId) {
        // overwrite
        h2t[hitId] = trackId;
      }
    }

    seedIdx++;
  }

  for (std::uint32_t trackIdx = 0; trackIdx < vSeedCandidates.size();
       ++trackIdx) {
    const std::uint32_t nTotal = vSeedCandidates[trackIdx].spacePoints.size();
    std::uint32_t nOther = 0;

    const std::uint32_t trackId = vTrackIds[trackIdx];

    for (const auto& h : vSeedCandidates[trackIdx].spacePoints) {
      const std::uint32_t hitId = h + 1;

      const std::uint32_t tid = h2t[hitId];

      // taken by a better candidate
      if (tid != trackId) {
        nOther++;
      }
    }

    if (nOther > m_cfg.hitShareThreshold * nTotal) {
      // reject
      vSeedCandidates[trackIdx].isClone = -1;
    }
  }
}

}  // namespace Acts::Experimental
