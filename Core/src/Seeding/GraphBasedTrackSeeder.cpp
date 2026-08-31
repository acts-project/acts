// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"

#include "Acts/Seeding/GbtsTrackingFilter.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <numbers>
#include <stdexcept>
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
  // buildTheGraph pre-computes the loosest tau ratio threshold it can apply,
  // which assumes the correction only ever widens the cut.
  if (m_cfg.tauRatioCorr < 0) {
    throw std::invalid_argument(
        "GraphBasedTrackSeeder: tauRatioCorr must not be negative");
  }

  m_tauLut = parseTauLookupTable(m_cfg.lutInputFile);
}

GbtsNodeStorage GraphBasedTrackSeeder::makeNodeStorage() const {
  GbtsNodeStorage::Config config;
  config.useClusterWidthCuts = m_cfg.useClusterWidthCuts;
  config.maxEndcapClusterWidth = m_cfg.maxEndcapClusterWidth;
  config.moduleHalfLengthY = m_cfg.moduleHalfLengthY;
  config.moduleEdgeTolerance = m_cfg.moduleEdgeTolerance;
  config.phiSliceWidth = m_cfg.phiSliceWidth;

  return GbtsNodeStorage(config, m_geometry, m_tauLut);
}

void GraphBasedTrackSeeder::createSeeds(const SpacePointContainer& spacePoints,
                                        const GbtsRoiDescriptor& roi,
                                        const GbtsTrackingFilter& filter,
                                        const Options& options,
                                        SeedContainer& outputSeeds) const {
  GbtsNodeStorage nodeStorage = makeNodeStorage();

  const auto layerColumn = spacePoints.column<std::uint32_t>("layerId");
  const auto clusterWidthColumn = spacePoints.column<float>("clusterWidth");
  const auto localPositionColumn = spacePoints.column<float>("localPositionY");

  nodeStorage.extend(spacePoints, layerColumn, clusterWidthColumn,
                     localPositionColumn);

  nodeStorage.finalize();

  createSeeds(nodeStorage, roi, filter, options, outputSeeds);
}

void GraphBasedTrackSeeder::createSeeds(GbtsNodeStorage& nodeStorage,
                                        const GbtsRoiDescriptor& roi,
                                        const GbtsTrackingFilter& filter,
                                        const Options& options,
                                        SeedContainer& outputSeeds) const {
  ACTS_DEBUG("Loaded " << nodeStorage.numberOfNodes() << " graph nodes");

  std::vector<detail::GbtsEdge> edgeStorage;

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
  extractSeedsFromTheGraph(maxLevel, graphStats.first, nodeStorage, edgeStorage,
                           vOutputSeeds, filter);

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

detail::GbtsTauLookupTable GraphBasedTrackSeeder::parseTauLookupTable(
    const std::string& lutInputFile) {
  if (!m_cfg.useClusterWidthCuts) {
    return {};
  }
  if (lutInputFile.empty()) {
    throw std::runtime_error("Cannot find tau lookup table file");
  }

  std::ifstream ifs(std::string(lutInputFile).c_str());
  if (!ifs.is_open()) {
    throw std::runtime_error("Failed to open tau lookup table file");
  }

  detail::GbtsTauLookupTable tauLut;
  tauLut.reserve(100);

  // per line: cluster width, bulk tau bounds, near-edge tau bounds. The width
  // is dropped - rows are located by index, never searched.
  float clusterWidth{};
  detail::GbtsTauBounds bounds;
  while (ifs >> clusterWidth >> bounds.minTau >> bounds.maxTau >>
         bounds.minTauNearEdge >> bounds.maxTauNearEdge) {
    tauLut.push_back(bounds);
  }

  if (!ifs.eof()) {
    // ended if parse error present, not clean EOF
    throw std::runtime_error("Stopped reading LUT file due to parse error");
  }

  return tauLut;
}

std::pair<std::int32_t, std::int32_t> GraphBasedTrackSeeder::buildTheGraph(
    const GbtsRoiDescriptor& roi, GbtsNodeStorage& nodeStorage,
    std::vector<detail::GbtsEdge>& edgeStorage, const Options& options) const {
  // used to calculate Z cut on doublets
  const float cutZMinU =
      m_cfg.minZ0 + m_cfg.maxOuterRadius * static_cast<float>(roi.dzdrMin());
  const float cutZMaxU =
      m_cfg.maxZ0 + m_cfg.maxOuterRadius * static_cast<float>(roi.dzdrMax());

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

  // the loosest tau ratio threshold the triplet matching can apply
  const float maxTauRatioCut =
      m_cfg.tauRatioCut + (m_cfg.useAdaptiveCuts ? m_cfg.tauRatioCorr : 0.0f) +
      (nodeStorage.hasStrips() ? m_cfg.tauRatioCorrStrip : 0.0f);

  // the default sliding window along phi
  const float deltaPhi0 = 0.5f * m_cfg.phiSliceWidth;

  std::uint32_t nConnections = 0;

  edgeStorage.reserve(m_cfg.nMaxEdges);

  std::uint32_t nEdges = 0;

  // scale factor to get indexes of binned beamspot
  // assuming 16-bit z0 bitmask

  const std::uint32_t zBins = 16;
  const float z0HistoCoeff = zBins / (m_cfg.maxZ0 - m_cfg.minZ0 + 1e-6);

  const detail::GbtsNodeView nodeView = nodeStorage.nodeView();
  const std::span<const detail::GbtsNodeParams> params =
      nodeStorage.nodeParams();
  const std::span<detail::GbtsNodeEdgeInfo> edgeInfo =
      nodeStorage.nodeEdgeInfo();

  // reused across bin groups so that the windows are allocated once
  std::vector<SlidingWindow> phiSlidingWindow;

  const bool calibrate = m_cfg.calibrateStrips && nodeStorage.hasStrips();

  // Put a strip node back where a direction says it crossed. Both ends, since
  // the nominal position sits on the inner strip and the calibrated one on the
  // outer. A direction that misses either strip is no crossing at all.
  const auto calibrateNode = [&](const SpacePointIndex node,
                                 const std::array<float, 3>& direction,
                                 float& r, float& z) {
    std::array<float, 3> point{};
    if (!Acts::detail::calibrateOuterStripSpacePoint(
            direction, nodeStorage.strip(node), point,
            m_cfg.stripLengthTolerance)) {
      return false;
    }
    r = fastHypot(point[0], point[1]);
    z = point[2];
    return true;
  };

  // loop over bin groups
  for (const auto& bg : m_geometry->binGroups()) {
    const detail::GbtsEtaBinInfo& B1 = nodeStorage.etaBin(bg.first);

    if (B1.empty()) {
      continue;
    }

    const float rb1 = B1.minRadius;

    const std::uint32_t layerId1 = B1.layerId;

    const bool isBarrel1 = B1.type == GbtsLayerType::Barrel;
    const bool isPixel1 = B1.technology == GbtsLayerTechnology::Pixel;

    // prepare a sliding window for each non-empty bin2 in the group

    phiSlidingWindow.clear();

    // loop over n2 eta-bins in L2 layers
    for (const auto& b2Idx : bg.second) {
      const detail::GbtsEtaBinInfo& B2 = nodeStorage.etaBin(b2Idx);

      if (B2.empty()) {
        continue;
      }

      const float rb2 = B2.maxRadius;

      float deltaPhi = deltaPhi0;  // the default

      // override the default window width
      if (m_cfg.useEtaBinning) {
        const float absDr = std::fabs(rb2 - rb1);
        if (m_cfg.useOldTunings) {
          deltaPhi = m_cfg.minDeltaPhi + dPhiCoeff * absDr;
        } else {
          if (absDr < 60.0) {
            deltaPhi = 0.002f + 4.33e-4f * ptScale * absDr;
          } else {
            deltaPhi = 0.015f + 2.2e-4f * ptScale * absDr;
          }
        }
      }

      SlidingWindow& window = phiSlidingWindow.emplace_back();
      window.phiNodes = B2.phiNodes.data();
      window.numPhiNodes = static_cast<std::uint32_t>(B2.phiNodes.size());
      window.deltaPhi = deltaPhi;
      window.layerId = B2.layerId;
      window.type = B2.type;
      window.technology = B2.technology;
    }

    // in GBTSv3 the outer loop goes over n1 nodes in the Layer 1 bin
    for (SpacePointIndex n1Idx = B1.nodes.first; n1Idx < B1.nodes.second;
         ++n1Idx) {
      // initialization using the top watermark of the edge storage
      edgeInfo[n1Idx].firstEdge = nEdges;

      // the counter for the incoming graph edges created for n1
      std::uint16_t numCreatedEdges = 0;

      bool isConnected = false;

      std::array<std::uint8_t, 16> z0Histo = {};

      const detail::GbtsNodeParams& n1pars = params[n1Idx];

      const float phi1 = n1pars.phi;
      const float r1 = n1pars.r;
      const float z1 = n1pars.z;

      // the chord of a pair, which only the strip path below reads
      float x1 = 0.f;
      float y1 = 0.f;
      if (calibrate) {
        const std::array<float, 4>& position = nodeView.positions[n1Idx];
        x1 = position[0];
        y1 = position[1];
      }

      // the intermediate loop over sliding windows
      for (auto& slw : phiSlidingWindow) {
        const std::uint32_t lk2 = slw.layerId;

        const bool isBarrel2 = slw.type == GbtsLayerType::Barrel;

        const bool isPixel2 = slw.technology == GbtsLayerTechnology::Pixel;
        const bool stripPair = calibrate && (!isPixel1 || !isPixel2);

        const float deltaPhi = slw.deltaPhi;

        // sliding window phi1 +/- deltaPhi

        const float minPhi = phi1 - deltaPhi;
        const float maxPhi = phi1 + deltaPhi;

        // the inner loop over n2 nodes using sliding window
        for (std::uint32_t n2PhiIdx = slw.firstIt; n2PhiIdx < slw.numPhiNodes;
             ++n2PhiIdx) {
          const float phi2 = slw.phiNodes[n2PhiIdx].first;

          if (phi2 < minPhi) {
            // update the window position
            slw.firstIt = n2PhiIdx;
            continue;
          }
          if (phi2 > maxPhi) {
            // break and go to the next window
            break;
          }

          const SpacePointIndex n2Idx = slw.phiNodes[n2PhiIdx].second;

          const detail::GbtsNodeEdgeInfo& n2Info = edgeInfo[n2Idx];

          const std::uint16_t nodeInfo = n2Info.isConnected;

          // skip isolated nodes as their incoming edges lead to nowhere
          if ((layerId1 == 80000) && (nodeInfo == 0)) {
            continue;
          }

          const std::uint32_t n2FirstEdge = n2Info.firstEdge;
          const std::uint16_t n2NumEdges = n2Info.numEdges;
          const std::uint32_t n2LastEdge = n2FirstEdge + n2NumEdges;

          const detail::GbtsNodeParams& n2pars = params[n2Idx];

          const float r2 = n2pars.r;

          float dr = r2 - r1;

          // On the nominal radii, so an endcap pair the slide would have
          // opened up is lost here. The cheap reject is worth more.
          if (dr < m_cfg.minDeltaRadius) {
            continue;
          }

          const float z2 = n2pars.z;

          // the ends as the pair puts them: nominal, or slid along the strip
          // when resolved. Azimuth is kept as it was -- exactly right in the
          // barrel, and to the stereo angle in the endcap, where the strip
          // slid along is a hair off radial.
          float r1c = r1;
          float z1c = z1;
          float r2c = r2;
          float z2c = z2;

          if (stripPair) {
            // A circle's chord bisects the tangent directions, so turning the
            // chord by the pair's own azimuth difference gives each tangent.
            const std::array<float, 4>& position = nodeView.positions[n2Idx];
            const float dx = position[0] - x1;
            const float dy = position[1] - y1;
            const float dzc = z2 - z1;
            const float turn = phi2 - phi1;

            bool resolved = true;
            if (!isPixel1) {
              resolved = calibrateNode(
                  n1Idx, {dx + dy * turn, dy - dx * turn, dzc}, r1c, z1c);
            }
            if (resolved && !isPixel2) {
              resolved = calibrateNode(
                  n2Idx, {dx - dy * turn, dy + dx * turn, dzc}, r2c, z2c);
            }
            if (!resolved) {
              continue;
            }

            dr = r2c - r1c;
            if (dr < m_cfg.minDeltaRadius) {
              continue;
            }
          }

          const float dz = z2c - z1c;
          const float tau = dz / dr;
          const float ftau = std::fabs(tau);
          if (ftau > m_cfg.maxAbsTau) {
            continue;
          }

          if (ftau < n1pars.minTau) {
            continue;
          }
          if (ftau > n1pars.maxTau) {
            continue;
          }

          if (ftau < n2pars.minTau) {
            continue;
          }
          if (ftau > n2pars.maxTau) {
            continue;
          }

          const float z0 = z1c - r1c * tau;

          if (layerId1 == 80000) {  // check against non-empty z0 histogram
            if (!checkZ0BitMask(nodeInfo, z0, m_cfg.minZ0, z0HistoCoeff)) {
              continue;
            }
          }

          if (m_cfg.doubletFilterRZ) {
            if (z0 < m_cfg.minZ0 || z0 > m_cfg.maxZ0) {
              continue;
            }

            const float zOuter = z0 + m_cfg.maxOuterRadius * tau;

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

          const float hypotTau = fastHypot(1, tau);
          const float expEta = hypotTau - tau;
          // 1 / expEta, since (hypotTau - tau) * (hypotTau + tau) == 1. The
          // sum is also the better conditioned form for large tau.
          const float invExpEta = hypotTau + tau;

          // match edge candidate against edges incoming to n2
          if (m_cfg.matchBeforeCreate &&
              (layerId1 == 80000 || layerId1 == 81000)) {
            // we must have enough incoming edges to decide
            bool isGood = n2NumEdges <= 2;

            if (!isGood) {
              const float uat1 = invExpEta;

              for (std::uint32_t n2InIdx = n2FirstEdge; n2InIdx < n2LastEdge;
                   ++n2InIdx) {
                const float tau2 = edgeStorage[n2InIdx].p[0];
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

          const float dPhi2 = curv * r2c;
          const float dPhi1 = curv * r1c;

          if (nEdges < m_cfg.nMaxEdges) {
            edgeStorage.emplace_back(n1Idx, n2Idx, lk2, slw.type, expEta, curv,
                                     phi1 + dPhi1);

            ++numCreatedEdges;

            const std::uint32_t outEdgeIdx = nEdges;

            const float uat2 = invExpEta;
            const float phi2u = phi2 + dPhi2;
            const float curv2 = curv;

            // looking for neighbours of the new edge
            for (std::uint32_t inEdgeIdx = n2FirstEdge; inEdgeIdx < n2LastEdge;
                 ++inEdgeIdx) {
              detail::GbtsEdge* pS = &edgeStorage[inEdgeIdx];

              const float absTauRatio = std::abs(pS->p[0] * uat2 - 1.0f);

              // rejects most candidates before the layer bookkeeping below
              if (absTauRatio > maxTauRatioCut) {
                continue;
              }

              if (pS->nNei >= detail::kGbtsMaxEdgeNeighbours) {
                continue;
              }

              const std::uint32_t lk3 = pS->n2LayerId;

              const bool isBarrel3 = pS->n2Type == GbtsLayerType::Barrel;

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
              // The two doublets sharing a strip node resolved it separately,
              // so a triplet through a strip may disagree on tau by more. Any
              // of the three: the outer two carry their end's error into tau.
              if (m_cfg.tauRatioCorrStrip > 0.f &&
                  (!isPixel1 || !isPixel2 ||
                   nodeView.strip(pS->n2) != nullptr)) {
                addTauRatioCorr += m_cfg.tauRatioCorrStrip;
              }

              // bad match
              if (absTauRatio > m_cfg.tauRatioCut + addTauRatioCorr) {
                continue;
              }

              float dPhi = phi2u - pS->p[2];

              if (dPhi < -std::numbers::pi_v<float>) {
                dPhi += 2 * std::numbers::pi_v<float>;
              } else if (dPhi > std::numbers::pi_v<float>) {
                dPhi -= 2 * std::numbers::pi_v<float>;
              }

              if (std::abs(dPhi) > m_cfg.cutDPhiMax) {
                continue;
              }

              const float dcurv = curv2 - pS->p[1];

              if (dcurv < -m_cfg.cutDCurvMax || dcurv > m_cfg.cutDCurvMax) {
                continue;
              }

              // final check: cuts on pT and d0
              if (m_cfg.validateTriplets) {
                // Pixel barrel
                if (isBarrel1 && isBarrel2 && isBarrel3) {
                  const std::array<SpacePointIndex, 3> candidateTriplet = {
                      n1Idx, n2Idx, pS->n2};

                  if (!validateTriplet(nodeView, candidateTriplet, tripletPtMin,
                                       absTauRatio, m_cfg.tauRatioCut,
                                       options)) {
                    continue;
                  }
                }
              }

              pS->vNei[pS->nNei] = outEdgeIdx;
              ++pS->nNei;

              isConnected = true;  // there is at least one good match

              // edge confirmed - update z0 histogram

              const std::uint32_t z0BinIndex =
                  static_cast<std::uint32_t>(z0HistoCoeff * (z0 - m_cfg.minZ0));

              ++z0Histo[z0BinIndex];

              nConnections++;
            }
            nEdges++;
          }
        }  // loop over n2 (outer) nodes inside a sliding window on n2 bin
      }  // loop over sliding windows associated with n2 bins

      // updating the n1 node attributes

      edgeInfo[n1Idx].numEdges = numCreatedEdges;
      if (isConnected) {
        std::uint16_t z0BitMask = 0x0;

        for (std::uint32_t bIdx = 0; bIdx < 16; ++bIdx) {
          if (z0Histo[bIdx] == 0) {
            continue;
          }

          z0BitMask |= (1 << bIdx);
        }

        // non-zero mask indicates that there is at least one connected edge
        edgeInfo[n1Idx].isConnected = z0BitMask;
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
    const std::uint32_t nEdges,
    std::vector<detail::GbtsEdge>& edgeStorage) const {
  constexpr std::uint32_t maxIter = 15;

  std::int32_t maxLevel = 0;

  std::uint32_t iter = 0;

  std::vector<detail::GbtsEdge*> vOld;

  for (std::uint32_t edgeIndex = 0; edgeIndex < nEdges; ++edgeIndex) {
    detail::GbtsEdge* pS = &(edgeStorage[edgeIndex]);
    if (pS->nNei == 0) {
      continue;
    }

    // TODO: increment level for segments as they already have at least one
    // neighbour
    vOld.push_back(pS);
  }

  std::vector<detail::GbtsEdge*> vNew;
  vNew.reserve(vOld.size());

  // generate proposals
  for (; iter < maxIter; iter++) {
    vNew.clear();

    for (detail::GbtsEdge* pS : vOld) {
      std::int32_t nextLevel = pS->level;

      for (std::uint32_t nIdx = 0; nIdx < pS->nNei; ++nIdx) {
        const std::uint32_t nextEdgeIdx = pS->vNei[nIdx];

        detail::GbtsEdge* pN = &(edgeStorage[nextEdgeIdx]);

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
    std::uint32_t maxLevel, std::uint32_t nEdges,
    const GbtsNodeStorage& nodeStorage,
    std::vector<detail::GbtsEdge>& edgeStorage,
    std::vector<OutputSeedProperties>& vOutputSeeds,
    const GbtsTrackingFilter& filter) const {
  const detail::GbtsNodeView nodeView = nodeStorage.nodeView();
  // a triplet + 1 confirmation
  std::uint8_t minLevel = 3;

  if (m_cfg.lrtMode) {
    // a triplet + no confirmation
    minLevel = 2;
  }

  if (maxLevel < minLevel) {
    return;
  }

  std::vector<detail::GbtsEdge*> vChainHeads;

  vChainHeads.reserve(nEdges / 2);

  for (std::uint32_t edgeIndex = 0; edgeIndex < nEdges; ++edgeIndex) {
    detail::GbtsEdge* pS = &edgeStorage[edgeIndex];

    if (m_cfg.lrtMode || !m_cfg.addTriplets) {
      if (pS->level < minLevel) {
        continue;
      }
    } else {  // eta-dependent cut
      const float edgeAbsEta = std::abs(-std::log(pS->p[0]));

      if (edgeAbsEta > m_cfg.maxAbsEtaAddTripelts) {
        if (pS->level < minLevel) {
          continue;
        }
      } else {
        if (pS->level < minLevel - 1) {
          continue;
        }
      }
    }

    vChainHeads.push_back(pS);
  }

  if (vChainHeads.empty()) {
    return;
  }

  std::ranges::sort(vChainHeads, std::ranges::greater{},
                    [](const detail::GbtsEdge* e) { return e->level; });

  // backtracking

  std::vector<SeedCandidateProperties> vSeedCandidates;

  vSeedCandidates.reserve(vChainHeads.size());

  std::vector<std::pair<float, std::uint32_t>> vArgSort;

  vArgSort.reserve(vChainHeads.size());

  std::uint32_t seedCounter = 0;

  GbtsTrackingFilter::State filterState{};

  for (detail::GbtsEdge* pS : vChainHeads) {
    if (pS->level == -1) {
      continue;
    }

    detail::GbtsEdgeState rs =
        filter.followTrack(filterState, nodeView, edgeStorage, *pS);

    if (!rs.initialized) {
      continue;
    }

    const float seedAbsEta = std::abs(-std::log(pS->p[0]));

    const std::uint32_t chainLength = static_cast<std::uint32_t>(rs.vs.size());

    if (m_cfg.lrtMode || !m_cfg.addTriplets) {
      if (chainLength < minLevel) {
        continue;
      }
    } else {
      if (seedAbsEta > m_cfg.maxAbsEtaAddTripelts) {
        if (chainLength < minLevel) {
          continue;
        }
      } else {
        if (chainLength < static_cast<std::uint32_t>(minLevel) - 1u) {
          continue;
        }
      }
    }

    std::vector<SpacePointIndex> vN;

    for (auto sIt = rs.vs.rbegin(); sIt != rs.vs.rend(); ++sIt) {
      if (seedAbsEta > m_cfg.edgeMaskMinEta) {
        // mark as collected
        (*sIt)->level = -1;
      }

      if (sIt == rs.vs.rbegin()) {
        vN.push_back((*sIt)->n1);
      }

      vN.push_back((*sIt)->n2);
    }

    // a triplet is accepted if it makes it up to this point
    if (vN.size() < 3) {
      continue;
    }

    const std::uint32_t origSeedSize = vN.size();

    const float origSeedQuality = -rs.j / origSeedSize;

    std::uint32_t seedSplitFlag = (seedAbsEta < m_cfg.maxSeedSplitEta) &&
                                          (origSeedSize > 3) &&
                                          (origSeedSize <= 5)
                                      ? 1
                                      : 0;

    // split the seed by dropping spacepoints
    if (seedSplitFlag != 0) {
      // 2. "drop-outs" and the original seed candidate
      std::array<std::array<SpacePointIndex, 3>, 3> triplets{};

      // triplet parameter estimate
      std::array<float, 3> invRads{};

      triplets[0] = {vN[0], vN[origSeedSize / 2], vN[origSeedSize - 1]};

      // all but the first one
      const std::vector<SpacePointIndex> dropOut1(vN.begin() + 1, vN.end());

      triplets[1] = {dropOut1[0], dropOut1[(origSeedSize - 1) / 2],
                     dropOut1[origSeedSize - 2]};

      std::vector<SpacePointIndex> dopOut2;

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
        invRads[k] = estimateCurvature(nodeView, triplets[k]);
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

  // hit to track associations, indexed by graph node index
  std::vector<std::uint32_t> h2t(nodeStorage.numberOfNodes() + 1, 0);

  std::uint32_t trackId = 0;

  for (const auto& args : vArgSort) {
    const auto& seed = vSeedCandidates[args.second];
    ++trackId;

    // loop over space points indices
    for (const SpacePointIndex node : seed.nodes) {
      const std::uint32_t hitId = node + 1;

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
    const auto& seed = vSeedCandidates[args.second].nodes;

    const std::uint32_t nTotal = seed.size();

    std::uint32_t nOther = 0;

    trackId = trackIdx + 1;

    ++trackIdx;

    for (const SpacePointIndex node : seed) {
      const std::uint32_t hitId = node + 1;

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

    const auto& vN = seed.nodes;

    if (seed.forSeedSplitting == 0) {
      // add seed to output

      std::vector<std::uint32_t> vSpIdx;

      vSpIdx.resize(vN.size());

      for (std::uint32_t k = 0; k < vSpIdx.size(); k++) {
        vSpIdx[k] = nodeStorage.spacePointIndex(vN[k]);
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

        newSeed.emplace_back(nodeStorage.spacePointIndex(vN[k]));
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
    const detail::GbtsNodeView& nodeView,
    const std::array<SpacePointIndex, 3>& nodes) const {
  // conformal mapping with the center at the last spacepoint

  std::array<float, 2> u{};
  std::array<float, 2> v{};

  const detail::GbtsNodeProxy n0 = nodeView[nodes[2]];

  const float x0 = n0.x();
  const float y0 = n0.y();

  const float r0 = n0.r();

  const float cosA = x0 / r0;

  const float sinA = y0 / r0;

  for (std::uint32_t k = 0; k < 2; k++) {
    const detail::GbtsNodeProxy nk = nodeView[nodes[k]];

    const float dx = nk.x() - x0;

    const float dy = nk.y() - y0;

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
    const detail::GbtsNodeView& nodeView,
    const std::array<SpacePointIndex, 3>& candidateTriplet,
    const float tripletMinPt, const float tauRatio, const float tauRatioCut,
    const Options& options) const {
  // conformal mapping with the center at the middle spacepoint

  std::array<float, 2> u{};
  std::array<float, 2> v{};

  const detail::GbtsNodeProxy n0 = nodeView[candidateTriplet[1]];

  const float x0 = n0.x();
  const float y0 = n0.y();

  const float r0 = n0.r();

  const float cosA = x0 / r0;

  const float sinA = y0 / r0;

  for (std::uint32_t k = 0; k < 2; k++) {
    const std::uint32_t spIdx = (k == 1) ? 2 : k;

    const detail::GbtsNodeProxy nk = nodeView[candidateTriplet[spIdx]];

    const float dx = nk.x() - x0;

    const float dy = nk.y() - y0;

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
