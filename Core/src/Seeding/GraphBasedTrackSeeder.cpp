// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"

#include "Acts/Seeding/GbtsGraphBuilder.hpp"
#include "Acts/Seeding/GbtsTrackingFilter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Acts::Experimental {

GraphBasedTrackSeeder::DerivedConfig::DerivedConfig(const Config& config)
    : Config(config) {
  phiSliceWidth = 2 * std::numbers::pi_v<float> / config.nMaxPhiSlice;
}

GraphBasedTrackSeeder::GraphBasedTrackSeeder(
    const DerivedConfig& config, std::shared_ptr<GbtsGeometry> geometry,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(config),
      m_geometry(std::move(geometry)),
      m_logger(std::move(logger)) {
  if (m_cfg.phiSortBuckets > GbtsNodeStorage::kMaxPhiSortBuckets) {
    throw std::invalid_argument(
        "GraphBasedTrackSeeder: phiSortBuckets exceeds the maximum");
  }

  if (m_cfg.useClusterWidthCuts && m_cfg.tauLookupTable.empty()) {
    throw std::invalid_argument(
        "GraphBasedTrackSeeder: the cluster width cuts need a tau lookup "
        "table");
  }
}

GbtsNodeStorage GraphBasedTrackSeeder::makeNodeStorage() const {
  GbtsNodeStorage::Config config;
  config.useClusterWidthCuts = m_cfg.useClusterWidthCuts;
  config.maxEndcapClusterWidth = m_cfg.maxEndcapClusterWidth;
  config.moduleHalfLengthY = m_cfg.moduleHalfLengthY;
  config.moduleEdgeTolerance = m_cfg.moduleEdgeTolerance;
  config.phiSliceWidth = m_cfg.phiSliceWidth;
  config.phiIndexMargin = m_cfg.phiIndexMargin;
  config.phiSortBuckets = m_cfg.phiSortBuckets;
  config.tauLutBinWidth = m_cfg.tauLutBinWidth;

  return GbtsNodeStorage(config, m_geometry, m_cfg.tauLookupTable);
}

void GraphBasedTrackSeeder::createSeeds(const SpacePointContainer& spacePoints,
                                        const GbtsRoiDescriptor& roi,
                                        const GbtsGraphBuilder& graphBuilder,
                                        const GbtsTrackingFilter& filter,
                                        const Options& options,
                                        SeedContainer& outputSeeds) const {
  GbtsNodeStorage nodeStorage = makeNodeStorage();

  const auto layerColumn = spacePoints.column<GbtsLayerIndex>("gbtsLayerIndex");
  const auto clusterWidthColumn = spacePoints.column<float>("clusterWidth");
  const auto localPositionColumn = spacePoints.column<float>("localPositionY");

  nodeStorage.extend(spacePoints, layerColumn, clusterWidthColumn,
                     localPositionColumn);

  nodeStorage.finalize();

  createSeeds(nodeStorage, roi, graphBuilder, filter, options, outputSeeds);
}

void GraphBasedTrackSeeder::createSeeds(GbtsNodeStorage& nodeStorage,
                                        const GbtsRoiDescriptor& roi,
                                        const GbtsGraphBuilder& graphBuilder,
                                        const GbtsTrackingFilter& filter,
                                        const Options& options,
                                        SeedContainer& outputSeeds) const {
  ACTS_DEBUG("Loaded " << nodeStorage.numberOfNodes() << " graph nodes");

  GbtsGraph graph =
      graphBuilder.buildTheGraph(roi, nodeStorage, options.bFieldInZ);

  ACTS_DEBUG("Created graph with " << graph.nEdges << " edges and "
                                   << graph.nConnections << " edge links");

  if (graph.nEdges == 0 || graph.nConnections == 0) {
    ACTS_WARNING("Missing edges or edge connections");
  }

  const std::uint32_t maxLevel = graphBuilder.runCCA(graph);

  const auto minLevel =
      static_cast<std::uint8_t>(graphBuilder.config().minSeedLevel);
  if (maxLevel < minLevel) {
    return;
  }

  ACTS_DEBUG("Reached Level " << maxLevel << " after GNN iterations");

  std::vector<detail::GbtsEdge*> vChainHeads =
      graphBuilder.extractChainHeads(graph);

  if (vChainHeads.empty()) {
    ACTS_WARNING("No chains passed minimum edge requirement");
    return;
  }

  std::vector<OutputSeedProperties> vOutputSeeds;
  extractSeedsFromTheGraph(nodeStorage, graph.edgeStorage, vOutputSeeds, filter,
                           vChainHeads, graphBuilder);

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

// TODO: fix the extractSeedsFromGraph function
void GraphBasedTrackSeeder::extractSeedsFromTheGraph(
    const GbtsNodeStorage& nodeStorage,
    std::vector<detail::GbtsEdge>& edgeStorage,
    std::vector<OutputSeedProperties>& vOutputSeeds,
    const GbtsTrackingFilter& filter,
    std::vector<detail::GbtsEdge*>& vChainHeads,
    const GbtsGraphBuilder& graphBuilder) const {
  const detail::GbtsNodeView nodeView = nodeStorage.nodeView();
  // the chain selection is the graph builder's, so that the chains it handed
  // back and the candidates built from them are cut the same way
  const GbtsGraphBuilder::Config& graphCfg = graphBuilder.config();
  const auto minLevel = static_cast<std::uint8_t>(graphCfg.minSeedLevel);
  // `addTriplets` accepts a chain one level short. Signed: an uncollected
  // edge sits at level -1 and `minSeedLevel` may be configured to 0.
  const int minLevelAddTriplets = int{minLevel} - 1;

  //===== everything to hear should be in the new class

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

    if (!graphCfg.addTriplets) {
      if (chainLength < minLevel) {
        continue;
      }
    } else {
      if (seedAbsEta > graphCfg.maxAbsEtaAddTriplets) {
        if (chainLength < minLevel) {
          continue;
        }
      } else {
        if (minLevelAddTriplets > 0 &&
            chainLength < static_cast<std::uint32_t>(minLevelAddTriplets)) {
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

    const auto origSeedSize = static_cast<std::uint32_t>(vN.size());

    const float origSeedQuality = -rs.j / origSeedSize;

    bool seedSplitFlag = (seedAbsEta < m_cfg.maxSeedSplitEta) &&
                         (origSeedSize >= m_cfg.minSplitSeedSize) &&
                         (origSeedSize <= m_cfg.maxSplitSeedSize);

    // split the seed by dropping spacepoints
    if (seedSplitFlag) {
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
        seedSplitFlag = false;  // reset the flag
      }
    }

    vSeedCandidates.emplace_back(origSeedQuality, false, vN, seedSplitFlag);

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

    const auto nTotal = static_cast<std::uint32_t>(seed.size());

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
      vSeedCandidates[args.second].isClone = true;  // reject
    }
  }
  vOutputSeeds.reserve(vSeedCandidates.size());

  // drop the clones and split seeds if need be

  for (const auto& args : vArgSort) {
    const auto& seed = vSeedCandidates[args.second];

    if (seed.isClone) {
      continue;  // identified as a clone of a better candidate
    }

    const auto& vN = seed.nodes;

    if (!seed.forSeedSplitting) {
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

    const auto seedSize = static_cast<std::uint32_t>(vN.size());

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

// this stays where it is
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

}  // namespace Acts::Experimental
