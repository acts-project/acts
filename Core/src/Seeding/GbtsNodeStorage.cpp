// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsNodeStorage.hpp"

#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/SpacePointFormation/detail/StripSpacePointCalibrationImpl.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <utility>

namespace Acts::Experimental {

GbtsNodeStorage::GbtsNodeStorage(const Config& config,
                                 std::shared_ptr<const GbtsGeometry> geometry,
                                 detail::GbtsTauLookupTable tauLut)
    : m_cfg(config),
      m_geometry(std::move(geometry)),
      m_tauLut(std::move(tauLut)),
      m_nodes(SpacePointColumns::CopiedFromIndex |
              SpacePointColumns::PackedXYZR) {
  m_etaBins.resize(m_geometry->numBins());
  m_stagedPerBin.resize(m_geometry->numBins());
}

std::optional<std::uint32_t> GbtsNodeStorage::insert(
    const SpacePointIndex index, const float x, const float y, const float z,
    const GbtsLayerIndex layerIndex, const float clusterWidth,
    const float localPositionY) {
  const float r = fastHypot(x, y);
  const float phi = std::atan2(y, x);
  return insert(index, x, y, z, r, phi, layerIndex, clusterWidth,
                localPositionY);
}

std::optional<std::uint32_t> GbtsNodeStorage::insert(
    const SpacePointIndex index, const float x, const float y, const float z,
    const float r, const float phi, const GbtsLayerIndex layerIndex,
    const float clusterWidth, const float localPositionY,
    const OuterStripSpacePointCalibrationDetails* strip) {
  const detail::GbtsLayer& layer = m_geometry->layerByIndex(layerIndex);
  const GbtsLayerDescription& description = layer.layerDescription();

  // wide pixel endcap clusters are dropped when the width cuts are on
  if (m_cfg.useClusterWidthCuts && description.type == GbtsLayerType::Endcap &&
      description.technology == GbtsLayerTechnology::Pixel &&
      clusterWidth > m_cfg.maxEndcapClusterWidth) {
    return std::nullopt;
  }

  const std::uint32_t bin = layer.getEtaBin(z, r);

  std::uint32_t stripIndex = detail::kNoStrip;
  // A strip pair on a pixel layer is never read: the seeder takes the strip
  // path per bin.
  if (strip != nullptr &&
      description.technology == GbtsLayerTechnology::Strip) {
    stripIndex = static_cast<std::uint32_t>(m_strips.size());
    // Derived once here rather than once per pair in the graph: it is six
    // cross products and a node takes part in many pairs.
    m_strips.push_back(
        Acts::detail::deriveOuterStripSpacePointCalibrationDetails(*strip));
  }

  m_stagedPerBin.at(bin).push_back(static_cast<std::uint32_t>(m_staged.size()));
  m_staged.emplace_back(index, x, y, z, r, phi, clusterWidth, localPositionY,
                        layerIndex, stripIndex);

  return bin;
}

void GbtsNodeStorage::extend(
    const SpacePointContainer& spacePoints,
    const ConstSpacePointColumnProxy<GbtsLayerIndex>& layerColumn,
    const ConstSpacePointColumnProxy<float>& clusterWidthColumn,
    const ConstSpacePointColumnProxy<float>& localPositionYColumn) {
  const bool strips =
      spacePoints.hasColumns(SpacePointColumns::StripCalibrationDetails);
  m_staged.reserve(m_staged.size() + spacePoints.size());
  for (const auto& sp : spacePoints) {
    insert(sp, layerColumn, clusterWidthColumn, localPositionYColumn, strips);
  }
}

std::vector<std::uint32_t> GbtsNodeStorage::sortBinByPhi(
    const std::vector<std::uint32_t>& staged) const {
  const std::uint32_t nBuckets = m_cfg.phiSortBuckets;
  std::array<std::vector<std::pair<float, std::uint32_t>>,
             kMaxPhiSortBuckets + 1>
      phiBuckets;

  for (const std::uint32_t stagedIdx : staged) {
    const float phi = m_staged[stagedIdx].phi;
    const auto bIdx = static_cast<std::uint32_t>(
        0.5 * nBuckets * (phi / std::numbers::pi_v<float> + 1.0f));
    phiBuckets[bIdx].emplace_back(phi, stagedIdx);
  }

  // Nodes with identical phi are ordered by insertion index.
  for (std::uint32_t bucket = 0; bucket <= nBuckets; ++bucket) {
    std::ranges::sort(phiBuckets[bucket]);
  }

  std::vector<std::uint32_t> sorted;
  sorted.reserve(staged.size());
  for (std::uint32_t bucket = 0; bucket <= nBuckets; ++bucket) {
    for (const auto& [phi, stagedIdx] : phiBuckets[bucket]) {
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
    detail::GbtsEtaBinInfo& binInfo = m_etaBins[bin];

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
    // every node in a bin is on the same layer, so any of them will do
    const GbtsLayerDescription& description =
        m_geometry->layerDescription(m_staged[staged.front()].layer);
    binInfo.layerId = description.id;
    binInfo.type = description.type;
    binInfo.technology = description.technology;
  }

  // Created now that the container has its final size, so that each column is
  // allocated in a single resize.
  m_paramsColumn.emplace(
      m_nodes.createColumn<detail::GbtsNodeParams>("gbtsNodeParams"));
  m_edgeInfoColumn.emplace(
      m_nodes.createColumn<detail::GbtsNodeEdgeInfo>("gbtsNodeEdgeInfo"));

  // Reordered into node order, so that the pairs a bin reads sit together;
  // left empty when nothing carries one, which is what `hasStrips` reports.
  std::vector<OuterStripSpacePointCalibrationDetailsDerived> strips;
  if (!m_strips.empty()) {
    m_stripIndex.assign(nodeOrder.size(), detail::kNoStrip);
    strips.reserve(m_strips.size());
  }

  std::span<detail::GbtsNodeParams> params = m_paramsColumn->data();
  for (std::uint32_t node = 0; node < nodeOrder.size(); ++node) {
    const StagedNode& staged = m_staged[nodeOrder[node]];

    if (staged.strip != detail::kNoStrip) {
      m_stripIndex[node] = static_cast<std::uint32_t>(strips.size());
      strips.push_back(m_strips[staged.strip]);
    }

    params[node].phi = staged.phi;
    params[node].r = staged.r;
    params[node].z = staged.z;

    // minTau and maxTau keep their "do not cut" defaults unless the lookup
    // table narrows them
    if (m_cfg.useClusterWidthCuts) {
      applyTauCuts(staged, params[node]);
    }
  }

  m_strips = std::move(strips);

  generatePhiIndexing(m_cfg.phiIndexMargin * m_cfg.phiSliceWidth);

  m_staged.clear();
  m_staged.shrink_to_fit();
  m_stagedPerBin.clear();
  m_stagedPerBin.shrink_to_fit();
}

void GbtsNodeStorage::applyTauCuts(const StagedNode& staged,
                                   detail::GbtsNodeParams& params) const {
  const GbtsLayerDescription& description =
      m_geometry->layerDescription(staged.layer);

  // the table is trained on pixel barrel clusters
  if (description.technology != GbtsLayerTechnology::Pixel ||
      description.type != GbtsLayerType::Barrel) {
    return;
  }

  // by the reciprocal, not the division: 1/0.05f is exactly 20, the division
  // is not, and the difference lands on the bin edges
  const auto lutBinIdx =
      static_cast<std::int32_t>(
          std::floor(staged.clusterWidth * (1.0f / m_cfg.tauLutBinWidth))) -
      1;

  if (lutBinIdx < 0 ||
      lutBinIdx >= static_cast<std::int32_t>(m_tauLut.size())) {
    return;
  }

  const detail::GbtsTauBounds& bounds =
      m_tauLut[static_cast<std::size_t>(lutBinIdx)];

  // close to the edge the cluster may be shortened, which the lookup table
  // covers with a separate pair of bounds
  const float dist2border =
      m_cfg.moduleHalfLengthY - std::abs(staged.localPositionY);
  const bool nearEdge = dist2border <= m_cfg.moduleEdgeTolerance;

  params.minTau = nearEdge ? bounds.minTauNearEdge : bounds.minTau;
  params.maxTau = nearEdge ? bounds.maxTauNearEdge : bounds.maxTau;

  if (params.maxTau < 0) {
    // insufficient training data, do not cut on tau
    params.maxTau = std::numeric_limits<float>::infinity();
  }
}

void GbtsNodeStorage::generatePhiIndexing(const float dphi) {
  const std::span<const detail::GbtsNodeParams> params = m_paramsColumn->data();

  for (detail::GbtsEtaBinInfo& bin : m_etaBins) {
    if (bin.empty()) {
      continue;
    }

    const SpacePointIndex begin = bin.nodes.first;
    const SpacePointIndex end = bin.nodes.second;

    for (SpacePointIndex node = begin; node < end; ++node) {
      const float phi = params[node].phi;
      if (phi <= std::numbers::pi_v<float> - dphi) {
        continue;
      }
      bin.phiNodes.emplace_back(phi - 2 * std::numbers::pi_v<float>, node);
    }

    for (SpacePointIndex node = begin; node < end; ++node) {
      bin.phiNodes.emplace_back(params[node].phi, node);
    }

    for (SpacePointIndex node = begin; node < end; ++node) {
      const float phi = params[node].phi;
      if (phi >= -std::numbers::pi_v<float> + dphi) {
        break;
      }
      bin.phiNodes.emplace_back(phi + 2 * std::numbers::pi_v<float>, node);
    }
  }
}

}  // namespace Acts::Experimental
