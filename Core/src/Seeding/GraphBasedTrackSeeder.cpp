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
#include <limits>
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

void GraphBasedTrackSeeder::createSeeds(
    GbtsNodeStorage& nodeStorage, const GbtsRoiDescriptor& roi,
    const GbtsGraphBuilder& graphBuilder, const GbtsTrackingFilter& filter,
    const Options& options, SeedContainer& outputSeeds,
    std::vector<std::optional<GbtsSeedFit>>* outputFits) const {
  ACTS_DEBUG("Loaded " << nodeStorage.numberOfNodes() << " graph nodes");

  GbtsGraph graph =
      graphBuilder.buildTheGraph(roi, nodeStorage, options.bFieldInZ);

  ACTS_DEBUG("Created graph with " << graph.edgeStorage.size() << " edges and "
                                   << graph.nConnections << " edge links");

  if (graph.edgeStorage.empty() || graph.nConnections == 0) {
    ACTS_WARNING("Missing edges or edge connections");
  }

  // the best path per edge settles its own levels
  if (!m_cfg.bestPathPerEdge) {
    const std::uint32_t maxLevel = graph.runCCA(m_cfg.ccaMaxIterations);

    ACTS_DEBUG("Reached Level " << maxLevel << " after GNN iterations");
  }

  std::vector<OutputSeedProperties> vOutputSeeds;
  extractSeedsFromTheGraph(nodeStorage, graph, vOutputSeeds, filter,
                           outputFits != nullptr);

  ACTS_DEBUG("GBTS created " << vOutputSeeds.size() << " seeds");
  if (vOutputSeeds.empty()) {
    ACTS_WARNING("No Seed Candidates");
  }

  // add to output seed container
  for (const auto& seed : vOutputSeeds) {
    auto newSeed = outputSeeds.createSeed();
    newSeed.assignSpacePointIndices(seed.spacePoints);
    newSeed.quality() = seed.seedQuality;
    if (outputFits != nullptr) {
      outputFits->push_back(seed.fit);
    }
  }
}

void GraphBasedTrackSeeder::extractSeedsFromTheGraph(
    const GbtsNodeStorage& nodeStorage, GbtsGraph& graph,
    std::vector<OutputSeedProperties>& vOutputSeeds,
    const GbtsTrackingFilter& filter, const bool withFits) const {
  const detail::GbtsNodeView nodeView = nodeStorage.nodeView();
  // `addTriplets` accepts a chain one level short. Signed: an uncollected
  // edge sits at level -1 and `minSeedLevel` may be configured to 0.
  const std::int8_t minLevelAddTriplets =
      static_cast<std::int8_t>(m_cfg.minSeedLevel - 1);

  std::vector<SeedCandidateProperties> vSeedCandidates;

  std::vector<std::pair<float, std::uint32_t>> vArgSort;

  std::uint32_t seedCounter = 0;

  // Turn the fitted chain of a head edge into a seed candidate, marking its
  // edges as collected if asked to.
  const auto addCandidate = [&](const detail::GbtsEdgeState& rs,
                                const detail::GbtsEdge& head,
                                const bool maskEdges) {
    const float seedAbsEta = std::abs(-std::log(head.p[0]));

    const std::uint32_t chainLength = static_cast<std::uint32_t>(rs.vs.size());

    if (!m_cfg.addTriplets) {
      if (chainLength < m_cfg.minSeedLevel) {
        return;
      }
    } else {
      if (seedAbsEta > m_cfg.maxAbsEtaAddTriplets) {
        if (chainLength < m_cfg.minSeedLevel) {
          return;
        }
      } else {
        if (minLevelAddTriplets > 0 &&
            chainLength < static_cast<std::uint32_t>(minLevelAddTriplets)) {
          return;
        }
      }
    }

    std::vector<SpacePointIndex> vN;

    for (auto sIt = rs.vs.rbegin(); sIt != rs.vs.rend(); ++sIt) {
      if (maskEdges && seedAbsEta > m_cfg.edgeMaskMinEta) {
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
      return;
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

    if (withFits) {
      // filter state at the innermost node, x' along the first edge (inwards)
      const double slope2 = 1. + static_cast<double>(rs.x[1]) * rs.x[1];
      vSeedCandidates.back().fit = GbtsSeedFit{
          .position = Vector3(rs.refX * rs.c - rs.x[0] * rs.s,
                              rs.refX * rs.s + rs.x[0] * rs.c, rs.y[0]),
          // outwards: along -x', with the slope of the state
          .phi = std::atan2(-rs.s - rs.c * rs.x[1], -rs.c + rs.s * rs.x[1]),
          .cotTheta = rs.y[1],
          .curvature = rs.x[2] / (slope2 * std::sqrt(slope2))};
    }

    vArgSort.emplace_back(origSeedQuality, seedCounter);

    ++seedCounter;
  };

  if (m_cfg.bestPathPerEdge) {
    extractBestPathPerEdge(nodeView, graph, filter, addCandidate);
  } else {
    std::vector<detail::GbtsEdge*> vChainHeads = graph.extractChainHeads(
        m_cfg.minSeedLevel, m_cfg.addTriplets, m_cfg.maxAbsEtaAddTriplets);

    if (vChainHeads.empty()) {
      ACTS_WARNING("No chains passed minimum edge requirement");
      return;
    }

    // backtracking

    vSeedCandidates.reserve(vChainHeads.size());

    vArgSort.reserve(vChainHeads.size());

    GbtsTrackingFilter::State filterState{};

    for (detail::GbtsEdge* pS : vChainHeads) {
      if (pS->level == -1) {
        continue;
      }

      detail::GbtsEdgeState rs =
          filter.followTrack(filterState, nodeView, graph.edgeStorage, *pS);

      if (!rs.initialized) {
        continue;
      }

      addCandidate(rs, *pS, true);
    }
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
      vOutputSeeds.back().fit = seed.fit;

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

template <typename add_candidate_t>
void GraphBasedTrackSeeder::extractBestPathPerEdge(
    const detail::GbtsNodeView& nodeView, GbtsGraph& graph,
    const GbtsTrackingFilter& filter, add_candidate_t& addCandidate) const {
  std::vector<detail::GbtsEdge>& edgeStorage = graph.edgeStorage;
  const auto nEdges = static_cast<std::uint32_t>(edgeStorage.size());
  const auto minLevel = static_cast<std::uint32_t>(m_cfg.minSeedLevel);

  // The level of an edge is the length of the longest chain of edges going
  // outwards from it. `vNei` holds the edges further in, which were all
  // created after the edge, so a single pass in edge order settles it. Longer
  // chains than `ccaMaxIterations` edges are unsettled and take no part.
  const auto maxLength = static_cast<std::uint32_t>(m_cfg.ccaMaxIterations);
  const std::uint32_t unsettled = maxLength + 1;
  std::vector<std::uint32_t> level(nEdges, 1);
  for (std::uint32_t e = 0; e < nEdges; ++e) {
    const detail::GbtsEdge& edge = edgeStorage[e];
    const std::uint32_t next = std::min(level[e] + 1, unsettled);
    for (std::uint32_t k = 0; k < edge.nNei; ++k) {
      std::uint32_t& innerLevel = level[edge.vNei[k]];
      innerLevel = std::max(innerLevel, next);
    }
  }

  // A root has no settled edge further in. Every path goes from an edge
  // inwards to a root, each step to an edge one level up, i.e. along the
  // longest chain outwards of the root.
  std::vector<char> isRoot(nEdges, 1);
  for (std::uint32_t e = 0; e < nEdges; ++e) {
    const detail::GbtsEdge& edge = edgeStorage[e];
    for (std::uint32_t k = 0; k < edge.nNei; ++k) {
      if (level[edge.vNei[k]] <= maxLength) {
        isRoot[e] = 0;
        break;
      }
    }
  }

  // The longest path from an edge to a root, zero if none, to skip the edges
  // and branches that cannot give a long enough path. The edges further in
  // come later, so a single pass backwards settles it.
  std::vector<std::uint32_t> rootReach(nEdges, 0);
  for (std::uint32_t e = nEdges; e-- > 0;) {
    if (level[e] > maxLength) {
      continue;
    }
    if (isRoot[e] != 0) {
      rootReach[e] = 1;
      continue;
    }
    const detail::GbtsEdge& edge = edgeStorage[e];
    for (std::uint32_t k = 0; k < edge.nNei; ++k) {
      const std::uint32_t inner = edge.vNei[k];
      if (level[inner] == level[e] + 1 && rootReach[inner] > 0) {
        rootReach[e] = std::max(rootReach[e], rootReach[inner] + 1);
      }
    }
  }

  // Depth first search from every edge inwards, fitting on the way; the best
  // complete path wins the edge.
  struct Frame {
    std::uint32_t edge{};
    std::uint32_t next{};
    detail::GbtsEdgeState state;
  };
  std::vector<Frame> stack(maxLength);
  std::vector<std::uint32_t> path(maxLength);
  std::vector<std::uint32_t> bestPath(maxLength);

  for (std::uint32_t head = 0; head < nEdges; ++head) {
    detail::GbtsEdge& headEdge = edgeStorage[head];

    // `addTriplets` accepts a path one edge short at low eta
    std::uint32_t minLength = minLevel;
    if (m_cfg.addTriplets && minLength > 0 &&
        std::abs(-std::log(headEdge.p[0])) <= m_cfg.maxAbsEtaAddTriplets) {
      --minLength;
    }

    if (level[head] > maxLength || rootReach[head] < minLength) {
      continue;
    }

    std::uint32_t depth = 0;
    Frame& first = stack[0];
    first.edge = head;
    first.next = 0;
    first.state.initialize(headEdge, nodeView, filter.m_cfg.initialVarianceX,
                           filter.m_cfg.initialVarianceY);
    if (!filter.update(nodeView, headEdge, first.state)) {
      continue;
    }

    float bestQuality = std::numeric_limits<float>::max();
    std::uint32_t bestLength = 0;
    detail::GbtsEdgeState bestState;

    while (true) {
      Frame& frame = stack[depth];
      const detail::GbtsEdge& edge = edgeStorage[frame.edge];

      if (frame.next == 0 && isRoot[frame.edge] != 0) {
        // a complete path, scored as the seed quality
        const std::uint32_t length = depth + 1;
        const float quality = -frame.state.j / static_cast<float>(length + 1);
        if (length >= minLength && quality < bestQuality) {
          bestQuality = quality;
          bestLength = length;
          bestState = frame.state;
          for (std::uint32_t d = 0; d <= depth; ++d) {
            bestPath[d] = stack[d].edge;
          }
        }
      }

      bool descended = false;
      while (frame.next < edge.nNei) {
        const std::uint32_t inner = edge.vNei[frame.next++];
        // one level up, towards a root far enough in
        if (level[inner] != level[frame.edge] + 1 || rootReach[inner] == 0 ||
            depth + 1 + rootReach[inner] < minLength) {
          continue;
        }
        Frame& nextFrame = stack[depth + 1];
        nextFrame.state = frame.state;
        if (!filter.update(nodeView, edgeStorage[inner], nextFrame.state)) {
          continue;
        }
        nextFrame.edge = inner;
        nextFrame.next = 0;
        ++depth;
        descended = true;
        break;
      }
      if (descended) {
        continue;
      }
      if (depth == 0) {
        break;
      }
      --depth;
    }

    if (bestLength == 0) {
      continue;
    }

    bestState.vs.clear();
    for (std::uint32_t d = 0; d < bestLength; ++d) {
      bestState.vs.push_back(&edgeStorage[bestPath[d]]);
    }
    addCandidate(bestState, headEdge, false);
  }
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

}  // namespace Acts::Experimental
