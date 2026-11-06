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
#include <span>
#include <tuple>
#include <utility>

#include "detail/SortByKey.hpp"

namespace Acts::Experimental {

namespace {

/// The phi sort's result and the storage `detail::sortByKey` expects
/// (`scratch`, `buckets`), kept from bin to bin.
struct PhiSortBuffers {
  /// A bin's (phi, index into its staged nodes) pairs in phi order
  std::vector<std::pair<float, std::uint32_t>> order;
  std::vector<std::pair<float, std::uint32_t>> scratch;
  std::vector<std::uint32_t> buckets;
};

/// Fill `buffers.order` with a bin's (phi, staged index) pairs in phi order.
template <typename staged_node_t>
void sortBinByPhi(std::span<const staged_node_t> staged,
                  PhiSortBuffers& buffers) {
  using Item = std::pair<float, std::uint32_t>;

  buffers.order.clear();
  for (std::uint32_t stagedIdx = 0; stagedIdx < staged.size(); ++stagedIdx) {
    buffers.order.emplace_back(staged[stagedIdx].phi, stagedIdx);
  }

  // Sort by phi. Nodes with the same phi are ordered by z, then r, so the
  // result does not depend on the order the space points were inserted in.
  // Nodes at the same position stay in insertion order, as the sort is
  // stable.
  const auto phiOf = [](const Item& item) { return item.first; };
  const auto byPosition = [&staged](const Item& a, const Item& b) {
    const staged_node_t& nodeA = staged[a.second];
    const staged_node_t& nodeB = staged[b.second];
    return std::tie(nodeA.z, nodeA.r) < std::tie(nodeB.z, nodeB.r);
  };

  Acts::detail::sortByKey(std::span(buffers.order), buffers, phiOf, byPosition);
}

/// Fill a bin's phi index: its nodes in phi order, those within `dphi` of the
/// wrap-around repeated shifted by 2pi, so the sliding window never wraps.
void fillPhiNodes(detail::GbtsEtaBinInfo& binInfo,
                  std::span<const std::pair<float, std::uint32_t>> phiOrder,
                  const float dphi) {
  constexpr float pi = std::numbers::pi_v<float>;
  constexpr float twoPi = 2 * pi;

  // nodes [0, head) lie within dphi above -pi, nodes [tail, n) within dphi
  // below pi (the bin is in phi order)
  const auto n = static_cast<std::uint32_t>(phiOrder.size());
  std::uint32_t head = 0;
  while (head < n && phiOrder[head].first < -pi + dphi) {
    ++head;
  }
  std::uint32_t tail = n;
  while (tail > 0 && phiOrder[tail - 1].first > pi - dphi) {
    --tail;
  }
  const SpacePointIndex first = binInfo.nodes.first;
  binInfo.phiNodes.reserve(n + head + (n - tail));
  for (std::uint32_t k = tail; k < n; ++k) {
    binInfo.phiNodes.emplace_back(phiOrder[k].first - twoPi, first + k);
  }
  for (std::uint32_t k = 0; k < n; ++k) {
    binInfo.phiNodes.emplace_back(phiOrder[k].first, first + k);
  }
  for (std::uint32_t k = 0; k < head; ++k) {
    binInfo.phiNodes.emplace_back(phiOrder[k].first + twoPi, first + k);
  }
}

}  // namespace

GbtsNodeStorage::GbtsNodeStorage(const Config& config,
                                 std::shared_ptr<const GbtsGeometry> geometry,
                                 GbtsTauLookupTable tauLut)
    : m_cfg(config),
      m_geometry(std::move(geometry)),
      m_tauLut(std::move(tauLut)),
      m_etaBins(m_geometry->numBins()),
      m_stagedPerBin(m_geometry->numBins()) {}

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

  m_stagedPerBin.at(bin).emplace_back(index, x, y, z, r, phi, clusterWidth,
                                      localPositionY, layerIndex, stripIndex);

  return bin;
}

void GbtsNodeStorage::extend(
    const SpacePointContainer& spacePoints,
    const ConstSpacePointColumnProxy<GbtsLayerIndex>& layerColumn,
    const ConstSpacePointColumnProxy<float>& clusterWidthColumn,
    const ConstSpacePointColumnProxy<float>& localPositionYColumn) {
  const bool strips =
      spacePoints.hasColumns(SpacePointColumns::StripCalibrationDetails);
  for (const auto& sp : spacePoints) {
    insert(sp, layerColumn, clusterWidthColumn, localPositionYColumn, strips);
  }
}

void GbtsNodeStorage::finalize() {
  std::uint32_t nNodes = 0;
  for (const std::vector<StagedNode>& staged : m_stagedPerBin) {
    nNodes += static_cast<std::uint32_t>(staged.size());
  }

  m_positions.reserve(nNodes);
  m_spacePointIndex.reserve(nNodes);
  m_layers.reserve(nNodes);
  m_params.reserve(nNodes);
  // filled while the graph is built
  m_edgeInfo.resize(nNodes);

  const float dphi = m_cfg.phiIndexMargin * m_cfg.phiSliceWidth;

  // Reordered into node order, so that the pairs a bin reads sit together;
  // left empty when nothing carries one, which is what `hasStrips` reports.
  std::vector<OuterStripSpacePointCalibrationDetailsDerived> strips;
  if (!m_strips.empty()) {
    m_stripIndex.assign(nNodes, detail::kNoStrip);
    strips.reserve(m_strips.size());
  }

  PhiSortBuffers phiSort;
  for (std::uint32_t bin = 0; bin < m_etaBins.size(); ++bin) {
    detail::GbtsEtaBinInfo& binInfo = m_etaBins[bin];

    binInfo.nodes = {numberOfNodes(), numberOfNodes()};

    const std::vector<StagedNode>& staged = m_stagedPerBin[bin];
    if (staged.empty()) {
      continue;
    }

    sortBinByPhi(std::span(staged), phiSort);
    const std::vector<std::pair<float, std::uint32_t>>& phiOrder =
        phiSort.order;

    for (const auto& [phi, stagedIdx] : phiOrder) {
      const StagedNode& node = staged[stagedIdx];
      const std::uint32_t nodeIdx = numberOfNodes();

      m_positions.push_back({node.x, node.y, node.z, node.r});
      m_spacePointIndex.push_back(node.spacePointIndex);
      m_layers.push_back(node.layer);

      if (node.strip != detail::kNoStrip) {
        m_stripIndex[nodeIdx] = static_cast<std::uint32_t>(strips.size());
        strips.push_back(m_strips[node.strip]);
      }

      // minTau and maxTau keep their "do not cut" defaults unless the lookup
      // table narrows them
      detail::GbtsNodeParams params{.phi = phi, .r = node.r, .z = node.z};
      if (m_cfg.useClusterWidthCuts) {
        applyTauCuts(node, params);
      }
      m_params.push_back(params);

      binInfo.minRadius = std::min(binInfo.minRadius, node.r);
      binInfo.maxRadius = std::max(binInfo.maxRadius, node.r);
      binInfo.minZ = std::min(binInfo.minZ, node.z);
      binInfo.maxZ = std::max(binInfo.maxZ, node.z);
    }

    binInfo.nodes.second = numberOfNodes();
    // every node in a bin is on the same layer, so any of them will do
    const GbtsLayerDescription& description =
        m_geometry->layerDescription(staged.front().layer);
    binInfo.barrelOrder = description.barrelOrder;
    binInfo.type = description.type;
    binInfo.technology = description.technology;

    fillPhiNodes(binInfo, phiOrder, dphi);
  }

  m_strips = std::move(strips);

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

  const GbtsTauBounds& bounds = m_tauLut[static_cast<std::size_t>(lutBinIdx)];

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

}  // namespace Acts::Experimental
