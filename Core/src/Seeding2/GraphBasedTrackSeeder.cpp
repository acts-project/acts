// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding2/GraphBasedTrackSeeder.hpp"

#include "Acts/Seeding2/GbtsTrackingFilter.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <numbers>
#include <utility>
#include <vector>

namespace Acts::Experimental {

GraphBasedTrackSeeder::DerivedConfig::DerivedConfig(const Config& config)
    : Config(config) {
  phiSliceWidth = 2 * std::numbers::pi_v<float> / config.nMaxPhiSlice;
}

GraphBasedTrackSeeder::Options::Options(float bFieldInZ_)
    : bFieldInZ(bFieldInZ_) {
  ptCoeff = 0.5f * bFieldInZ * Acts::UnitConstants::m;
}

GraphBasedTrackSeeder::GraphBasedTrackSeeder(
    const DerivedConfig& config, std::shared_ptr<GbtsGeometry> geometry,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(config),
      m_geometry(std::move(geometry)),
      m_logger(std::move(logger)) {
  m_mlLut = parseGbtsMlLookupTable(m_cfg.lutInputFile);
}

void GraphBasedTrackSeeder::createSeeds(const SpacePointContainer2& spacePoints,
                                        const GbtsRoiDescriptor& roi,
                                        const std::uint32_t maxLayers,
                                        const GbtsTrackingFilter& filter,
                                        const Options& options,
                                        SeedContainer2& outputSeeds) const {
  GbtsNodeStorage nodeStorage(m_geometry, m_mlLut);

  std::vector<std::vector<GbtsNode>> nodesPerLayer =
      createNodes(spacePoints, maxLayers);
  std::uint32_t nPixelLoaded = 0;
  std::uint32_t nStripLoaded = 0;

  for (std::uint16_t l = 0; l < nodesPerLayer.size(); l++) {
    const std::vector<GbtsNode>& nodes = nodesPerLayer[l];

    if (nodes.empty()) {
      continue;
    }

    const bool isPixel = true;
    // placeholder for now until strip hits are added in
    if (isPixel) {
      nPixelLoaded += nodeStorage.loadPixelGraphNodes(
          l, nodes, m_cfg.useMl, m_cfg.maxEndcapClusterWidth);
    } else {
      nStripLoaded += nodeStorage.loadStripGraphNodes(l, nodes);
    }
  }
  ACTS_DEBUG("Loaded " << nPixelLoaded << " pixel space points and "
                       << nStripLoaded << " strip space points");

  nodeStorage.sortByPhi();

  nodeStorage.initializeNodes(m_cfg.useMl);

  nodeStorage.generatePhiIndexing(1.5f * m_cfg.phiSliceWidth);

  std::vector<GbtsEdge> edgeStorage;

  std::pair<std::int32_t, std::int32_t> graphStats =
      buildTheGraph(roi, nodeStorage, edgeStorage, options);

  ACTS_DEBUG("Created graph with " << graphStats.first << " edges and "
                                   << graphStats.second << " edge links");

  if (graphStats.first == 0 || graphStats.second == 0) {
    ACTS_WARNING("Missing edges or edge connections");
  }

  std::uint32_t maxLevel = runCCA(graphStats.first, edgeStorage);

  ACTS_DEBUG("Reached Level " << maxLevel << " after GNN iterations");

  std::vector<OutputSeedProperties> vOutputSeeds;
  extractSeedsFromTheGraph(maxLevel, graphStats.first, spacePoints.size(),
                           edgeStorage, vOutputSeeds, filter);

  ACTS_DEBUG("GBTS created " << vOutputSeeds.size() << " seeds");
  if (vOutputSeeds.empty()) {
    ACTS_WARNING("No Seed Candidates");
  }

  // add to output seed container
  for (const auto& seed : vOutputSeeds) {
    auto newSeed = outputSeeds.createSeed();
    newSeed.assignSpacePointIndices(seed.spacePoints);
    newSeed.quality() = seed.seedQuality;
  }
}

GbtsMlLookupTable GraphBasedTrackSeeder::parseGbtsMlLookupTable(
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

std::vector<std::vector<GbtsNode>> GraphBasedTrackSeeder::createNodes(
    const SpacePointContainer2& spacePoints,
    const std::uint32_t maxLayers) const {
  auto layerColumn = spacePoints.column<std::uint32_t>("layerId");
  auto clusterWidthColumn = spacePoints.column<float>("clusterWidth");
  auto localPositionColumn = spacePoints.column<float>("localPositionY");

  std::vector<std::vector<GbtsNode>> nodesPerLayer(maxLayers);
  // reserve for better efficiency
  for (auto& v : nodesPerLayer) {
    v.reserve(10000);
  }

  for (const auto& sp : spacePoints) {
    // for every sp in container,
    // add its variables to nodeStorage organised by layer
    const std::uint16_t layer = sp.extra(layerColumn);

    // add node to storage
    GbtsNode& node = nodesPerLayer[layer].emplace_back(layer);

    // fill the node with space point variables

    node.x = sp.x();
    node.y = sp.y();
    node.z = sp.z();
    node.r = sp.r();
    node.phi = sp.phi();
    node.idx = sp.index();
    node.pcw = sp.extra(clusterWidthColumn);
    node.locPosY = sp.extra(localPositionColumn);
  }

  return nodesPerLayer;
}

std::pair<std::int32_t, std::int32_t> GraphBasedTrackSeeder::buildTheGraph(
    const GbtsRoiDescriptor& roi, GbtsNodeStorage& nodeStorage,
    std::vector<GbtsEdge>& edgeStorage, const Options& options) const {
  // phi cut for triplets
  const float cutDPhiMax = m_cfg.lrtMode ? 0.07f : 0.012f;
  // curv cut for triplets
  const float cutDCurvMax = m_cfg.lrtMode ? 0.015f : 0.001f;
  // tau cut for doublets and triplets
  const float cutTauRatioMax = m_cfg.lrtMode ? 0.015f : m_cfg.tauRatioCut;
  const float minZ0 = m_cfg.lrtMode ? -600.0f : static_cast<float>(roi.zMin());
  const float maxZ0 = m_cfg.lrtMode ? 600.0f : static_cast<float>(roi.zMax());
  const float minDeltaPhi = m_cfg.lrtMode ? 0.01f : 0.001f;

  // used to calculate Z cut on doublets
  const float maxOuterRadius = m_cfg.lrtMode ? 1050.0f : 550.0f;

  const float cutZMinU =
      minZ0 + maxOuterRadius * static_cast<float>(roi.dzdrMin());
  const float cutZMaxU =
      maxZ0 + maxOuterRadius * static_cast<float>(roi.dzdrMax());

  // correction due to limited pT resolution
  const float tripletPtMin = 0.8f * m_cfg.minPt;

  // to re-scale original tunings done for the 900 MeV pT cut
  const float ptScale = 0.9f / m_cfg.minPt;

  const float maxCurv = options.ptCoeff / tripletPtMin;

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
  const float deltaPhi0 = 0.5f * m_cfg.phiSliceWidth;

  std::uint32_t nConnections = 0;

  edgeStorage.reserve(m_cfg.nMaxEdges);

  std::uint32_t nEdges = 0;

  // scale factor to get indexes of binned beamspot
  // assuming 16-bit z0 bitmask

  const std::uint32_t zBins = 16;
  const float z0HistoCoeff = zBins / (maxZ0 - minZ0 + 1e-6);

  // loop over bin groups
  for (const auto& bg : m_geometry->binGroups()) {
    GbtsEtaBin& B1 = nodeStorage.getEtaBin(bg.first);

    if (B1.empty()) {
      continue;
    }

    const float rb1 = B1.minRadius;

    const std::uint32_t layerId1 = B1.layerId;

    const bool isBarrel1 = (layerId1 / 10000) == 8;

    // prepare a sliding window for each bin2 in the group

    std::vector<SlidingWindow> phiSlidingWindow;

    // initialization using default ctor
    phiSlidingWindow.resize(bg.second.size());

    std::uint32_t winIdx = 0;

    // loop over n2 eta-bins in L2 layers
    for (const auto& b2Idx : bg.second) {
      const GbtsEtaBin& B2 = nodeStorage.getEtaBin(b2Idx);

      if (B2.empty()) {
        ++winIdx;
        continue;
      }

      const float rb2 = B2.maxRadius;

      float deltaPhi = deltaPhi0;  // the default

      // override the default window width
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

      phiSlidingWindow[winIdx].bin = &B2;
      phiSlidingWindow[winIdx].hasNodes = true;
      phiSlidingWindow[winIdx].deltaPhi = deltaPhi;
      ++winIdx;
    }

    // in GBTSv3 the outer loop goes over n1 nodes in the Layer 1 bin
    for (std::uint32_t n1Idx = 0; n1Idx < B1.vn.size(); ++n1Idx) {
      // initialization using the top watermark of the edge storage
      B1.vFirstEdge[n1Idx] = nEdges;

      // the counter for the incoming graph edges created for n1
      std::uint16_t numCreatedEdges = 0;

      bool isConnected = false;

      std::array<std::uint8_t, 16> z0Histo = {};

      const std::array<float, 5>& n1pars = B1.params[n1Idx];

      const float phi1 = n1pars[2];
      const float r1 = n1pars[3];
      const float z1 = n1pars[4];

      // the intermediate loop over sliding windows
      for (auto& slw : phiSlidingWindow) {
        if (!slw.hasNodes) {
          continue;
        }

        const GbtsEtaBin& B2 = *slw.bin;

        const std::uint32_t lk2 = B2.layerId;

        const bool isBarrel2 = (lk2 / 10000) == 8;

        const float deltaPhi = slw.deltaPhi;

        // sliding window phi1 +/- deltaPhi

        const float minPhi = phi1 - deltaPhi;
        const float maxPhi = phi1 + deltaPhi;

        // the inner loop over n2 nodes using sliding window
        for (std::uint32_t n2PhiIdx = slw.firstIt;
             n2PhiIdx < B2.vPhiNodes.size(); ++n2PhiIdx) {
          const float phi2 = B2.vPhiNodes[n2PhiIdx].first;

          if (phi2 < minPhi) {
            // update the window position
            slw.firstIt = n2PhiIdx;
            continue;
          }
          if (phi2 > maxPhi) {
            // break and go to the next window
            break;
          }

          const std::uint32_t n2Idx = B2.vPhiNodes[n2PhiIdx].second;

          const std::uint16_t nodeInfo = B2.vIsConnected[n2Idx];

          // skip isolated nodes as their incoming edges lead to nowhere
          if ((layerId1 == 80000) && (nodeInfo == 0)) {
            continue;
          }

          const std::uint32_t n2FirstEdge = B2.vFirstEdge[n2Idx];
          const std::uint16_t n2NumEdges = B2.vNumEdges[n2Idx];
          const std::uint32_t n2LastEdge = n2FirstEdge + n2NumEdges;

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

          const float z0 = z1 - r1 * tau;

          if (layerId1 == 80000) {  // check against non-empty z0 histogram
            if (!checkZ0BitMask(nodeInfo, z0, minZ0, z0HistoCoeff)) {
              continue;
            }
          }

          if (m_cfg.doubletFilterRZ) {
            if (z0 < minZ0 || z0 > maxZ0) {
              continue;
            }

            const float zOuter = z0 + maxOuterRadius * tau;

            if (zOuter < cutZMinU || zOuter > cutZMaxU) {
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

          const float expEta = fastHypot(1, tau) - tau;

          // match edge candidate against edges incoming to n2
          if (m_cfg.matchBeforeCreate &&
              (layerId1 == 80000 || layerId1 == 81000)) {
            // we must have enough incoming edges to decide
            bool isGood = n2NumEdges <= 2;

            if (!isGood) {
              const float uat1 = 1.0f / expEta;

              for (std::uint32_t n2InIdx = n2FirstEdge; n2InIdx < n2LastEdge;
                   ++n2InIdx) {
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

            ++numCreatedEdges;

            const std::uint32_t outEdgeIdx = nEdges;

            const float uat2 = 1.f / expEta;
            const float phi2u = phi2 + dPhi2;
            const float curv2 = curv;

            // looking for neighbours of the new edge
            for (std::uint32_t inEdgeIdx = n2FirstEdge; inEdgeIdx < n2LastEdge;
                 ++inEdgeIdx) {
              GbtsEdge* pS = &(edgeStorage.at(inEdgeIdx));

              if (pS->nNei >= gbtsNumSegConns) {
                continue;
              }

              const std::uint32_t lk3 =
                  m_geometry->layerIdByIndex(pS->n2->layer);

              const bool isBarrel3 = (lk3 / 10000) == 8;

              const float absTauRatio = std::abs(pS->p[0] * uat2 - 1.0f);
              float addTauRatioCorr = 0;

              if (m_cfg.useAdaptiveCuts) {
                if (isBarrel1 && isBarrel2 && isBarrel3) {
                  const bool noGap =
                      ((lk3 - lk2) == 1000) && ((lk2 - layerId1) == 1000);

                  // assume more scattering due to the layer in between
                  if (!noGap) {
                    addTauRatioCorr = m_cfg.tauRatioCorr;
                  }
                } else {
                  bool mixedTriplet = isBarrel1 && isBarrel2 && !isBarrel3;
                  if (mixedTriplet) {
                    addTauRatioCorr = m_cfg.tauRatioCorr;
                  }
                }
              }
              // bad match
              if (absTauRatio > cutTauRatioMax + addTauRatioCorr) {
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

              // final check: cuts on pT and d0
              if (m_cfg.validateTriplets) {
                // Pixel barrel
                if (isBarrel1 && isBarrel2 && isBarrel3) {
                  const std::array<const GbtsNode*, 3> candidateTriplet = {
                      B1.vn[n1Idx], B2.vn[n2Idx], pS->n2};

                  if (!validateTriplet(candidateTriplet, tripletPtMin,
                                       absTauRatio, cutTauRatioMax, options)) {
                    continue;
                  }
                }
              }

              pS->vNei[pS->nNei] = outEdgeIdx;
              ++pS->nNei;

              isConnected = true;  // there is at least one good match

              // edge confirmed - update z0 histogram

              const std::uint32_t z0BinIndex =
                  static_cast<std::uint32_t>(z0HistoCoeff * (z0 - minZ0));

              ++z0Histo[z0BinIndex];

              nConnections++;
            }
            nEdges++;
          }
        }  // loop over n2 (outer) nodes inside a sliding window on n2 bin
      }  // loop over sliding windows associated with n2 bins

      // updating the n1 node attributes

      B1.vNumEdges[n1Idx] = numCreatedEdges;
      if (isConnected) {
        std::uint16_t z0BitMask = 0x0;

        for (std::uint32_t bIdx = 0; bIdx < 16; ++bIdx) {
          if (z0Histo[bIdx] == 0) {
            continue;
          }

          z0BitMask |= (1 << bIdx);
        }

        // non-zero mask indicates that there is at least one connected edge
        B1.vIsConnected[n1Idx] = z0BitMask;
      }

    }  // loop over n1 (inner) nodes
  }  // loop over bin groups: a single n1 bin and multiple n2 bins

  if (nEdges >= m_cfg.nMaxEdges) {
    ACTS_WARNING(
        "Maximum number of graph edges exceeded - possible efficiency loss "
        << nEdges);
  }
  return std::make_pair(nEdges, nConnections);
}

std::int32_t GraphBasedTrackSeeder::runCCA(
    const std::uint32_t nEdges, std::vector<GbtsEdge>& edgeStorage) const {
  constexpr std::uint32_t maxIter = 15;

  std::int32_t maxLevel = 0;

  std::uint32_t iter = 0;

  std::vector<GbtsEdge*> vOld;

  for (std::uint32_t edgeIndex = 0; edgeIndex < nEdges; ++edgeIndex) {
    GbtsEdge* pS = &(edgeStorage[edgeIndex]);
    if (pS->nNei == 0) {
      continue;
    }

    // TODO: increment level for segments as they already have at least one
    // neighbour
    vOld.push_back(pS);
  }

  std::vector<GbtsEdge*> vNew;
  vNew.reserve(vOld.size());

  // generate proposals
  for (; iter < maxIter; iter++) {
    vNew.clear();

    for (GbtsEdge* pS : vOld) {
      std::int32_t nextLevel = pS->level;

      for (std::uint32_t nIdx = 0; nIdx < pS->nNei; ++nIdx) {
        const std::uint32_t nextEdgeIdx = pS->vNei[nIdx];

        GbtsEdge* pN = &(edgeStorage[nextEdgeIdx]);

        if (pS->level == pN->level) {
          nextLevel = pS->level + 1;
          vNew.push_back(pS);
          break;
        }
      }

      // proposal
      pS->next = static_cast<std::int8_t>(nextLevel);
    }

    // update

    std::uint32_t nChanges = 0;

    for (auto pS : vNew) {
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

    vOld.swap(vNew);
    vNew.clear();
  }

  return maxLevel;
}

void GraphBasedTrackSeeder::extractSeedsFromTheGraph(
    std::uint32_t maxLevel, std::uint32_t nEdges, std::int32_t nHits,
    std::vector<GbtsEdge>& edgeStorage,
    std::vector<OutputSeedProperties>& vOutputSeeds,
    const GbtsTrackingFilter& filter) const {
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

  std::vector<SeedCandidateProperties> vSeedCandidates;

  vSeedCandidates.reserve(vSeeds.size());

  std::vector<std::pair<float, std::uint32_t>> vArgSort;

  vArgSort.reserve(vSeeds.size());

  std::uint32_t seedCounter = 0;

  GbtsTrackingFilter::State filterState{};

  for (GbtsEdge* pS : vSeeds) {
    if (pS->level == -1) {
      continue;
    }

    GbtsEdgeState rs = filter.followTrack(filterState, edgeStorage, *pS);

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

    const std::uint32_t origSeedSize = vN.size();

    const float origSeedQuality = -rs.j / origSeedSize;

    std::uint32_t seedSplitFlag = (seedEta < m_cfg.maxSeedSplitEta) &&
                                          (origSeedSize > 3) &&
                                          (origSeedSize <= 5)
                                      ? 1
                                      : 0;

    // split the seed by dropping spacepoints
    if (seedSplitFlag != 0) {
      // 2. "drop-outs" and the original seed candidate
      std::array<std::array<const GbtsNode*, 3>, 3> triplets{};

      // triplet parameter estimate
      std::array<float, 3> invRads{};

      triplets[0] = {vN[0], vN[origSeedSize / 2], vN[origSeedSize - 1]};

      // all but the first one
      const std::vector<const GbtsNode*> dropOut1(vN.begin() + 1, vN.end());

      triplets[1] = {dropOut1[0], dropOut1[(origSeedSize - 1) / 2],
                     dropOut1[origSeedSize - 2]};

      std::vector<const GbtsNode*> dopOut2;

      dopOut2.reserve(origSeedSize - 1);

      for (std::uint32_t k = 0; k < origSeedSize; k++) {
        if (k == origSeedSize / 2) {
          continue;  // drop the middle SP in the original seed
        }

        dopOut2.emplace_back(vN[k]);
      }

      triplets[2] = {dopOut2[0], dopOut2[(origSeedSize - 1) / 2],
                     dopOut2[origSeedSize - 2]};

      for (std::uint32_t k = 0; k < invRads.size(); k++) {
        invRads[k] = estimateCurvature(triplets[k]);
      }

      const std::array<float, 3> diffs = {std::abs(invRads[1] - invRads[0]),
                                          std::abs(invRads[2] - invRads[0]),
                                          std::abs(invRads[2] - invRads[1])};

      const bool confirmed = diffs[0] < m_cfg.maxInvRadDiff &&
                             diffs[1] < m_cfg.maxInvRadDiff &&
                             diffs[2] < m_cfg.maxInvRadDiff;

      if (confirmed) {
        seedSplitFlag = 0;  // reset the flag
      }
    }

    vSeedCandidates.emplace_back(origSeedQuality, 0, vN, seedSplitFlag);

    vArgSort.emplace_back(origSeedQuality, seedCounter);

    ++seedCounter;
  }

  // clone removal code goes below ...

  std::ranges::sort(vArgSort);

  std::vector<std::uint32_t> h2t(nHits + 1, 0);  // hit to track associations

  std::uint32_t trackId = 0;

  for (const auto& args : vArgSort) {
    const auto& seed = vSeedCandidates[args.second];
    ++trackId;

    // loop over space points indices
    for (const auto& h : seed.spacePoints) {
      const std::uint32_t hitId = h->idx + 1;

      const std::uint32_t tid = h2t[hitId];

      // unused hit or used by a lesser track
      if (tid == 0 || tid > trackId) {
        // overwrite
        h2t[hitId] = trackId;
      }
    }
  }

  std::uint32_t trackIdx = 0;

  for (const auto& args : vArgSort) {
    const auto& seed = vSeedCandidates[args.second].spacePoints;

    const std::uint32_t nTotal = seed.size();

    std::uint32_t nOther = 0;

    trackId = trackIdx + 1;

    ++trackIdx;

    for (const auto& h : seed) {
      const std::uint32_t hitId = h->idx + 1;

      const std::uint32_t tid = h2t[hitId];

      // taken by a better candidate
      if (tid != trackId) {
        nOther++;
      }
    }

    if (nOther > m_cfg.hitShareThreshold * nTotal) {
      // reject
      vSeedCandidates[args.second].isClone = -1;  // reject
    }
  }
  vOutputSeeds.reserve(vSeedCandidates.size());

  // drop the clones and split seeds if need be

  for (const auto& args : vArgSort) {
    const auto& seed = vSeedCandidates[args.second];

    if (seed.isClone != 0) {
      continue;  // identified as a clone of a better candidate
    }

    const auto& vN = seed.spacePoints;

    if (seed.forSeedSplitting == 0) {
      // add seed to output

      std::vector<std::uint32_t> vSpIdx;

      vSpIdx.resize(vN.size());

      for (std::uint32_t k = 0; k < vSpIdx.size(); k++) {
        vSpIdx[k] = vN[k]->idx;
      }

      vOutputSeeds.emplace_back(seed.seedQuality, vSpIdx);

      continue;
    }

    // seed split into "drop-out" seeds

    const std::uint32_t seedSize = vN.size();

    const std::array<std::size_t, 2> indices2drop = {
        0, seedSize / 2ul};  // the first and the middle

    for (const auto& skipIdx : indices2drop) {
      std::vector<std::uint32_t> newSeed;

      newSeed.reserve(seedSize - 1);

      for (std::uint32_t k = 0; k < seedSize; k++) {
        if (k == skipIdx) {
          continue;
        }

        newSeed.emplace_back(vN[k]->idx);
      }

      vOutputSeeds.emplace_back(seed.seedQuality, newSeed);
    }
  }
}

bool GraphBasedTrackSeeder::checkZ0BitMask(const std::uint16_t z0BitMask,
                                           const float z0, const float minZ0,
                                           const float z0HistoCoeff) const {
  if (z0BitMask == 0) {
    return true;
  }

  const float dz = z0 - minZ0;
  const std::int32_t z0BinIndex = static_cast<std::int32_t>(z0HistoCoeff * dz);

  if (((z0BitMask >> z0BinIndex) & 1) != 0) {
    return true;
  }

  // check adjacent bins as well

  const float z0Resolution = 2.5;

  const float dzm = dz - z0Resolution;

  std::int32_t nextBin = static_cast<std::int32_t>(z0HistoCoeff * dzm);

  if (nextBin >= 0 && nextBin != z0BinIndex) {
    if (((z0BitMask >> nextBin) & 1) != 0) {
      return true;
    }
  }

  const float dzp = dz + z0Resolution;

  nextBin = static_cast<std::int32_t>(z0HistoCoeff * dzp);

  if (nextBin < 16 && nextBin != z0BinIndex) {
    if (((z0BitMask >> nextBin) & 1) != 0) {
      return true;
    }
  }

  return false;
}

float GraphBasedTrackSeeder::estimateCurvature(
    const std::array<const GbtsNode*, 3>& nodes) const {
  // conformal mapping with the center at the last spacepoint

  std::array<float, 2> u{};
  std::array<float, 2> v{};

  const float x0 = nodes[2]->x;
  const float y0 = nodes[2]->y;

  const float r0 = nodes[2]->r;

  const float cosA = x0 / r0;

  const float sinA = y0 / r0;

  for (std::uint32_t k = 0; k < 2; k++) {
    const float dx = nodes[k]->x - x0;

    const float dy = nodes[k]->y - y0;

    const float r2Inv = 1.0 / (dx * dx + dy * dy);

    const float xn = dx * cosA + dy * sinA;

    const float yn = -dx * sinA + dy * cosA;

    u[k] = xn * r2Inv;
    v[k] = yn * r2Inv;
  }

  const float du = u[0] - u[1];

  if (du == 0.0) {
    return 0.0;
  }

  const float A = (v[0] - v[1]) / du;

  const float B = v[1] - A * u[1];

  // curavture in units of 1/mm
  return B / std::sqrt(1 + A * A);
}

bool GraphBasedTrackSeeder::validateTriplet(
    const std::array<const GbtsNode*, 3> candidateTriplet,
    const float tripletMinPt, const float tauRatio, const float tauRatioCut,
    const Options& options) const {
  // conformal mapping with the center at the middle spacepoint

  std::array<float, 2> u{};
  std::array<float, 2> v{};

  const float x0 = candidateTriplet[1]->x;
  const float y0 = candidateTriplet[1]->y;

  const float r0 = candidateTriplet[1]->r;

  const float cosA = x0 / r0;

  const float sinA = y0 / r0;

  for (std::uint32_t k = 0; k < 2; k++) {
    const std::uint32_t spIdx = (k == 1) ? 2 : k;

    const float dx = candidateTriplet[spIdx]->x - x0;

    const float dy = candidateTriplet[spIdx]->y - y0;

    const float r2Inv = 1.0f / (dx * dx + dy * dy);

    const float xn = dx * cosA + dy * sinA;

    const float yn = -dx * sinA + dy * cosA;

    u[k] = xn * r2Inv;
    v[k] = yn * r2Inv;
  }

  const float du = u[0] - u[1];

  if (du == 0.0) {
    return false;
  }

  const float A = (v[0] - v[1]) / du;

  const float B = v[1] - A * u[1];

  const float d0 = r0 * (B * r0 - A);

  if (std::abs(d0) > m_cfg.d0Max) {
    return false;
  }

  if (B != 0.0) {  // straight-line track is OK

    const float R = std::sqrt(1 + A * A) / B;

    // 1T magnetic field used
    const float pT = std::abs(options.bFieldInZ * R / 2);

    if (pT < tripletMinPt) {
      return false;
    }

    if (pT > 5 * tripletMinPt) {  // relatively high-pT track

      if (tauRatio > 0.9f * tauRatioCut) {
        return false;
      }
    }
  }

  return true;
}

}  // namespace Acts::Experimental
