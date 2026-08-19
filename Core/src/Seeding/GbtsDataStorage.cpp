// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsDataStorage.hpp"

#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

namespace Acts::Experimental {

GbtsNodeStorage::GbtsNodeStorage(Config config,
                                 std::shared_ptr<const GbtsGeometry> geometry,
                                 GbtsMlLookupTable mlLut)
    : m_cfg(std::move(config)),
      m_geometry(std::move(geometry)),
      m_mlLut(std::move(mlLut)),
      m_nodes(SpacePointColumns::CopiedFromIndex |
              SpacePointColumns::PackedXYZR) {
  m_etaBins.resize(m_geometry->numBins());
  m_stagedPerBin.resize(m_geometry->numBins());
}

std::optional<std::uint32_t> GbtsNodeStorage::insert(
    const SpacePointIndex index, const float x, const float y, const float z,
    const std::uint32_t layerIndex, const float clusterWidth,
    const float localPositionY) {
  const float r = fastHypot(x, y);
  const float phi = std::atan2(y, x);
  return insert(index, x, y, z, r, phi, layerIndex, clusterWidth,
                localPositionY);
}

std::optional<std::uint32_t> GbtsNodeStorage::insert(
    const SpacePointIndex index, const float x, const float y, const float z,
    const float r, const float phi, const std::uint32_t layerIndex,
    const float clusterWidth, const float localPositionY) {
  const GbtsLayer& layer = m_geometry->layerByIndex(layerIndex);

  const bool isBarrel = layer.layerDescription().type == GbtsLayerType::Barrel;

  // Pixel endcap nodes with a wide cluster are dropped when the machine
  // learning features are in use.
  if (m_cfg.useMl && !isBarrel && clusterWidth > m_cfg.maxEndcapClusterWidth &&
      layerIndex < m_cfg.isPixelLayer.size() &&
      m_cfg.isPixelLayer[layerIndex]) {
    return std::nullopt;
  }

  const std::int32_t binIndex = layer.getEtaBin(z, r);
  if (binIndex == -1) {
    return std::nullopt;
  }

  const auto bin = static_cast<std::uint32_t>(binIndex);

  m_stagedPerBin.at(bin).push_back(static_cast<std::uint32_t>(m_staged.size()));
  m_staged.push_back(StagedNode{index, x, y, z, r, phi, clusterWidth,
                                localPositionY,
                                static_cast<std::uint16_t>(layerIndex)});

  return bin;
}

void GbtsNodeStorage::extend(
    const SpacePointContainer& spacePoints,
    const ConstSpacePointColumnProxy<std::uint32_t>& layerColumn,
    const ConstSpacePointColumnProxy<float>& clusterWidthColumn,
    const ConstSpacePointColumnProxy<float>& localPositionYColumn) {
  m_staged.reserve(m_staged.size() + spacePoints.size());
  for (const auto& sp : spacePoints) {
    insert(sp, layerColumn, clusterWidthColumn, localPositionYColumn);
  }
}

std::vector<std::uint32_t> GbtsNodeStorage::sortBinByPhi(
    const std::vector<std::uint32_t>& staged) const {
  // TODO config
  constexpr std::uint32_t nBuckets = 31;
  std::array<std::vector<std::pair<float, std::uint32_t>>, 32> phiBuckets;

  for (const std::uint32_t stagedIdx : staged) {
    const float phi = m_staged[stagedIdx].phi;
    const auto bIdx = static_cast<std::uint32_t>(
        0.5 * nBuckets * (phi / std::numbers::pi_v<float> + 1.0f));
    phiBuckets[bIdx].emplace_back(phi, stagedIdx);
  }

  // Nodes with identical phi are ordered by insertion index.
  for (auto& bucket : phiBuckets) {
    std::ranges::sort(bucket);
  }

  std::vector<std::uint32_t> sorted;
  sorted.reserve(staged.size());
  for (const auto& bucket : phiBuckets) {
    for (const auto& [phi, stagedIdx] : bucket) {
      sorted.push_back(stagedIdx);
    }
  }

  return sorted;
}

void GbtsNodeStorage::finalize() {
  const auto nNodes = static_cast<std::uint32_t>(m_staged.size());

  m_nodes.reserve(nNodes);
  m_layers.reserve(nNodes);

  // Node order across all bins, used to fill the derived columns below.
  std::vector<std::uint32_t> nodeOrder;
  nodeOrder.reserve(nNodes);

  for (std::uint32_t bin = 0; bin < m_etaBins.size(); ++bin) {
    GbtsEtaBinInfo& binInfo = m_etaBins[bin];

    binInfo.nodes = {m_nodes.size(), m_nodes.size()};

    const std::vector<std::uint32_t>& staged = m_stagedPerBin[bin];
    if (staged.empty()) {
      continue;
    }

    const std::vector<std::uint32_t> sorted = sortBinByPhi(staged);

    float minRadius = std::numeric_limits<float>::max();
    float maxRadius = std::numeric_limits<float>::lowest();

    for (const std::uint32_t stagedIdx : sorted) {
      const StagedNode& node = m_staged[stagedIdx];

      auto newNode = m_nodes.createSpacePoint();
      newNode.copiedFromIndex() = node.spacePointIndex;
      newNode.xyzr() = std::array<float, 4>{node.x, node.y, node.z, node.r};

      m_layers.push_back(node.layer);
      nodeOrder.push_back(stagedIdx);

      minRadius = std::min(minRadius, node.r);
      maxRadius = std::max(maxRadius, node.r);
    }

    binInfo.nodes.second = m_nodes.size();
    binInfo.minRadius = minRadius;
    binInfo.maxRadius = maxRadius;
    binInfo.layerId =
        m_geometry->layerIdByIndex(m_staged[sorted.front()].layer);
  }

  // Created now that the container has its final size, so that each column is
  // allocated in a single resize.
  m_paramsColumn.emplace(
      m_nodes.createColumn<GbtsNodeParams>("gbtsNodeParams"));
  m_edgeInfoColumn.emplace(
      m_nodes.createColumn<GbtsNodeEdgeInfo>("gbtsNodeEdgeInfo"));

  std::span<GbtsNodeParams> params = m_paramsColumn->data();
  for (std::uint32_t node = 0; node < nodeOrder.size(); ++node) {
    const StagedNode& staged = m_staged[nodeOrder[node]];

    params[node].phi = staged.phi;
    params[node].r = staged.r;
    params[node].z = staged.z;

    // minTau and maxTau keep their "do not cut" defaults unless the lookup
    // table narrows them
    if (m_cfg.useMl) {
      applyMlTauCuts(staged, params[node]);
    }
  }

  generatePhiIndexing(1.5f * m_cfg.phiSliceWidth);

  m_staged.clear();
  m_staged.shrink_to_fit();
  m_stagedPerBin.clear();
  m_stagedPerBin.shrink_to_fit();
}

void GbtsNodeStorage::applyMlTauCuts(const StagedNode& staged,
                                     GbtsNodeParams& params) const {
  const GbtsLayer& layer = m_geometry->layerByIndex(staged.layer);

  // skip strips volumes: layers in range [1200X-1400X]
  if (layer.layerDescription().id < 20000) {
    return;
  }
  if (layer.layerDescription().type != GbtsLayerType::Barrel) {
    return;
  }

  // lut bin width is 0.05 mm
  const auto lutBinIdx =
      static_cast<std::int32_t>(std::floor(20 * staged.clusterWidth)) - 1;

  if (lutBinIdx < 0 || lutBinIdx >= static_cast<std::int32_t>(m_mlLut.size())) {
    return;
  }

  const GbtsMlLookupTable::value_type& lutBin = m_mlLut[lutBinIdx];

  // close to the edge the cluster may be shortened, which the lookup table
  // covers with a separate pair of bounds
  const float dist2border =
      m_cfg.moduleHalfLengthY - std::abs(staged.localPositionY);
  const bool nearEdge = dist2border <= m_cfg.moduleEdgeTolerance;

  params.minTau = nearEdge ? lutBin[3] : lutBin[1];
  params.maxTau = nearEdge ? lutBin[4] : lutBin[2];

  if (params.maxTau < 0) {
    // insufficient training data, do not cut on tau
    params.maxTau = std::numeric_limits<float>::infinity();
  }
}

void GbtsNodeStorage::generatePhiIndexing(const float dphi) {
  const std::span<const GbtsNodeParams> params = m_paramsColumn->data();

  for (GbtsEtaBinInfo& bin : m_etaBins) {
    if (bin.empty()) {
      continue;
    }

    const GbtsNodeIndex begin = bin.nodes.first;
    const GbtsNodeIndex end = bin.nodes.second;

    for (GbtsNodeIndex node = begin; node < end; ++node) {
      const float phi = params[node].phi;
      if (phi <= std::numbers::pi_v<float> - dphi) {
        continue;
      }
      bin.phiNodes.emplace_back(phi - 2 * std::numbers::pi_v<float>, node);
    }

    for (GbtsNodeIndex node = begin; node < end; ++node) {
      bin.phiNodes.emplace_back(params[node].phi, node);
    }

    for (GbtsNodeIndex node = begin; node < end; ++node) {
      const float phi = params[node].phi;
      if (phi >= -std::numbers::pi_v<float> + dphi) {
        break;
      }
      bin.phiNodes.emplace_back(phi + 2 * std::numbers::pi_v<float>, node);
    }
  }
}

}  // namespace Acts::Experimental
