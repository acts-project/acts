// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsGeometry.hpp"

#include "Acts/Utilities/MathHelpers.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <optional>
#include <unordered_map>

namespace Acts::Experimental::detail {

GbtsLayer::GbtsLayer(const GbtsLayerDescription& layerDescription,
                     const float etaBinWidth, const std::uint32_t bin0)
    : m_layerDescription(layerDescription) {
  float r1{};
  float r2{};
  float z1{};
  float z2{};
  if (m_layerDescription.type == GbtsLayerType::Barrel) {
    r1 = m_layerDescription.refCoord;
    r2 = m_layerDescription.refCoord;
    z1 = m_layerDescription.minBound;
    z2 = m_layerDescription.maxBound;
  } else if (m_layerDescription.type == GbtsLayerType::Endcap) {
    r1 = m_layerDescription.minBound;
    r2 = m_layerDescription.maxBound;
    z1 = m_layerDescription.refCoord;
    z2 = m_layerDescription.refCoord;
  } else {
    throw std::runtime_error("invalid layer type");
  }

  const float t1 = z1 / r1;
  const float eta1 = -std::log(fastHypot(1, t1) - t1);

  const float t2 = z2 / r2;
  const float eta2 = -std::log(fastHypot(1, t2) - t2);

  // increasing them slightly to avoid range_check exceptions
  const auto minEta = static_cast<float>(std::min(eta1, eta2) - 1e-6);
  const auto maxEta = static_cast<float>(std::max(eta1, eta2) + 1e-6);

  const float deltaEta = maxEta - minEta;

  // a layer thinner than one bin width keeps a single bin spanning it
  std::uint32_t numBins = 1;
  if (deltaEta >= etaBinWidth) {
    numBins = static_cast<std::uint32_t>(deltaEta / etaBinWidth);
    if (deltaEta - etaBinWidth * numBins > 0.5f * etaBinWidth) {
      ++numBins;
    }
  }

  m_binning = {bin0, numBins, minEta, deltaEta / static_cast<float>(numBins)};

  if (numBins == 1) {
    // the single bin spans the layer, so it takes the layer's own bounds
    if (m_layerDescription.type == GbtsLayerType::Barrel) {
      m_minRadius.push_back(m_layerDescription.refCoord - 2.0f);
      m_maxRadius.push_back(m_layerDescription.refCoord + 2.0f);
    } else {
      m_minRadius.push_back(m_layerDescription.minBound - 2.0f);
      m_maxRadius.push_back(m_layerDescription.maxBound + 2.0f);
    }
    m_minBinCoord.push_back(m_layerDescription.minBound);
    m_maxBinCoord.push_back(m_layerDescription.maxBound);
    return;
  }

  float eta = minEta + 0.5f * m_binning.etaBinWidth;

  for (std::uint32_t i = 0; i < numBins; ++i) {
    float e1 = eta - 0.5f * m_binning.etaBinWidth;
    float e2 = eta + 0.5f * m_binning.etaBinWidth;

    if (m_layerDescription.type == GbtsLayerType::Barrel) {
      m_minRadius.push_back(m_layerDescription.refCoord - 2.0f);
      m_maxRadius.push_back(m_layerDescription.refCoord + 2.0f);
      m_minBinCoord.push_back(m_layerDescription.refCoord * std::sinh(e1));
      m_maxBinCoord.push_back(m_layerDescription.refCoord * std::sinh(e2));
    } else {
      // for the positive endcap larger eta corresponds to smaller radius
      if (m_layerDescription.refCoord > 0) {
        std::swap(e1, e2);
      }
      float r = m_layerDescription.refCoord / std::sinh(e1);
      m_minBinCoord.push_back(r);
      m_minRadius.push_back(r - 2.0f);
      r = m_layerDescription.refCoord / std::sinh(e2);
      m_maxBinCoord.push_back(r);
      m_maxRadius.push_back(r + 2.0f);
    }

    eta += m_binning.etaBinWidth;
  }
}

bool GbtsLayer::checkCompatibility(const GbtsLayer& otherLayer,
                                   const std::uint32_t b1,
                                   const std::uint32_t b2, const float minZ0,
                                   const float maxZ0) const {
  const float z1min = m_minBinCoord.at(b1);
  const float z1max = m_maxBinCoord.at(b1);
  const float r1 = m_layerDescription.refCoord;

  const float tol = 5.0f;

  if (m_layerDescription.type == GbtsLayerType::Barrel &&
      otherLayer.m_layerDescription.type == GbtsLayerType::Barrel) {
    const float minB2 = otherLayer.m_minBinCoord.at(b2);
    const float maxB2 = otherLayer.m_maxBinCoord.at(b2);

    const float r2 = otherLayer.m_layerDescription.refCoord;

    const float A = r2 / (r2 - r1);
    const float B = r1 / (r2 - r1);

    const float z0Min = z1min * A - maxB2 * B;
    const float z0Max = z1max * A - minB2 * B;

    if (z0Max < minZ0 - tol || z0Min > maxZ0 + tol) {
      return false;
    }

    return true;
  }

  if (m_layerDescription.type == GbtsLayerType::Barrel &&
      otherLayer.m_layerDescription.type == GbtsLayerType::Endcap) {
    const float z2 = otherLayer.m_layerDescription.refCoord;
    const float r2max = otherLayer.m_maxBinCoord.at(b2);
    float r2min = otherLayer.m_minBinCoord.at(b2);

    if (r2max <= r1) {
      return false;
    }

    if (r2min <= r1) {
      r2min = r1 + 1e-3f;
    }

    float z0Max = 0;
    float z0Min = 0;

    if (z2 > 0) {
      z0Max = (z1max * r2max - z2 * r1) / (r2max - r1);
      z0Min = (z1min * r2min - z2 * r1) / (r2min - r1);
    } else {
      z0Max = (z1max * r2min - z2 * r1) / (r2min - r1);
      z0Min = (z1min * r2max - z2 * r1) / (r2max - r1);
    }

    if (z0Max < minZ0 - tol || z0Min > maxZ0 + tol) {
      return false;
    }
    return true;
  }

  if (m_layerDescription.type == GbtsLayerType::Endcap &&
      otherLayer.m_layerDescription.type == GbtsLayerType::Endcap) {
    const float z2 = otherLayer.m_layerDescription.refCoord;
    const float z1 = m_layerDescription.refCoord;
    const float r2max = otherLayer.m_maxBinCoord.at(b2);
    const float r2min = otherLayer.m_minBinCoord.at(b2);
    const float r1max = m_maxBinCoord.at(b1);
    const float r1min = m_minBinCoord.at(b1);

    if (r1min >= r2max) {
      return false;
    }

    if (z2 > 0) {  // positive endcap

      const float z0Max = z1 - r1min * (z2 - z1) / (r2max - r1min);

      if (z0Max < minZ0 - tol) {
        return false;
      }

      if (r2min > r1max) {
        const float z0Min = z1 - r1max * (z2 - z1) / (r2min - r1max);

        if (z0Min > maxZ0 + tol) {
          return false;
        }
      }
    } else {  // negative endcap
      const float z0Min = z1 - r1min * (z2 - z1) / (r2max - r1min);

      if (z0Min > maxZ0 + tol) {
        return false;
      }

      if (r2min > r1max) {
        const float z0Max = z1 - r1max * (z2 - z1) / (r2min - r1max);

        if (z0Max < minZ0 - tol) {
          return false;
        }
      }
    }
    return true;
  }

  if (m_layerDescription.type == GbtsLayerType::Endcap &&
      otherLayer.m_layerDescription.type == GbtsLayerType::Barrel) {
    const float z1 = m_layerDescription.refCoord;
    const float r1max = m_maxBinCoord.at(b1);
    const float r1min = m_minBinCoord.at(b1);

    const float z2min = otherLayer.m_minBinCoord.at(b2);
    const float z2max = otherLayer.m_maxBinCoord.at(b2);
    const float r2 = otherLayer.m_layerDescription.refCoord;

    if (r2 < r1min) {
      return false;
    }

    // interval 1

    float z0Min = z1 - (z2max - z1) / (r2 / r1max - 1);
    float z0Max = z1 - (z2max - z1) / (r2 / r1min - 1);

    if (z0Min > z0Max) {
      std::swap(z0Min, z0Max);
    }

    bool beyondRange = (z0Max < minZ0 - tol || z0Min > maxZ0 + tol);

    if (!beyondRange) {
      return true;
    }

    // interval 2

    z0Min = z1 - (z2min - z1) / (r2 / r1max - 1);
    z0Max = z1 - (z2min - z1) / (r2 / r1min - 1);

    if (z0Min > z0Max) {
      std::swap(z0Min, z0Max);
    }

    beyondRange = (z0Max < minZ0 - tol || z0Min > maxZ0 + tol);

    if (!beyondRange) {
      return true;
    }

    return false;
  }

  return true;
}

std::uint32_t GbtsLayer::getEtaBin(const float zh, const float rh) const {
  if (m_binning.numBins == 1) {
    return m_binning.firstBin;
  }

  const float t1 = zh / rh;
  const float eta = -std::log(fastHypot(1, t1) - t1);

  // signed: eta below the layer's range truncates negative, before the clamp
  std::int32_t idx = static_cast<std::int32_t>((eta - m_binning.minEta) /
                                               m_binning.etaBinWidth);
  if (idx < 0) {
    idx = 0;
  } else if (idx >= static_cast<std::int32_t>(m_binning.numBins)) {
    idx = static_cast<std::int32_t>(m_binning.numBins) - 1;
  }

  // index in the global storage
  return m_binning.firstBin + static_cast<std::uint32_t>(idx);
}

}  // namespace Acts::Experimental::detail

namespace Acts::Experimental {

namespace {
// key: bin. value: (outgoing, incoming) bins it connects to.
using BinConnections =
    std::unordered_map<std::uint32_t, std::pair<std::vector<std::uint32_t>,
                                                std::vector<std::uint32_t>>>;
}  // namespace

GbtsGeometry::GbtsGeometry(
    std::span<const GbtsLayerDescription> layerDescriptions,
    std::span<const GbtsLayerConnection> layerConnections,
    const float etaBinWidth, const GbtsZ0Range& z0Range, const Logger& logger)
    : m_etaBinWidth(etaBinWidth) {
  const float minZ0 = z0Range.min;
  const float maxZ0 = z0Range.max;

  for (const GbtsLayerDescription& layer : layerDescriptions) {
    const detail::GbtsLayer& pL = createLayer(layer, m_nEtaBins);
    m_nEtaBins += pL.binning().numBins;
  }

  // calculating bin tables in the connector...
  // calculate bin pairs for graph edge building

  std::optional<std::uint32_t> lastBin1;

  for (const GbtsLayerConnection& connection : layerConnections) {
    const detail::GbtsLayer* pL1 = layerById(connection.dst);  // n1
    const detail::GbtsLayer* pL2 = layerById(connection.src);  // n2
    if (pL1 == nullptr) {
      ACTS_WARNING("Skipping invalid dst layer " << connection.dst);
      continue;
    }
    if (pL2 == nullptr) {
      ACTS_WARNING("Skipping invalid src layer " << connection.src);
      continue;
    }

    const std::uint32_t nSrcBins = pL2->binning().numBins;
    const std::uint32_t nDstBins = pL1->binning().numBins;

    // loop over bins in Layer 1
    for (std::uint32_t b1 = 0; b1 < nDstBins; ++b1) {
      // loop over bins in Layer 2
      for (std::uint32_t b2 = 0; b2 < nSrcBins; ++b2) {
        if (!pL1->checkCompatibility(*pL2, b1, b2, minZ0, maxZ0)) {
          continue;
        }

        const std::uint32_t bin1Idx = pL1->binning().firstBin + b1;
        const std::uint32_t bin2Idx = pL2->binning().firstBin + b2;

        if (bin1Idx != lastBin1) {
          // adding a new group
          m_binGroups.push_back({bin1Idx, {bin2Idx}});
          lastBin1 = bin1Idx;
        } else {
          // extend the last group
          m_binGroups.back().links.push_back(bin2Idx);
        }
      }
    }
  }
  // find stages of eta-bin pairs using the graph ablation algorithm

  BinConnections binMap;

  // 1. create a map of bin-to-bin connections

  // initialize with empty links

  for (const auto& [bin1, bin2s] : m_binGroups) {
    // add to the map

    auto& bin1Links = binMap[bin1];

    for (const auto& bin2 : bin2s) {
      auto& bin2Links = binMap[bin2];

      bin1Links.second.push_back(bin2);  // incoming link bin1 <- bin2
      bin2Links.first.push_back(bin1);   // outgoing link bin2 -> bin1
    }
  }

  // copy bin map as original is still needed
  const BinConnections originalBinMap = binMap;

  // 2. find stages starting from the last one (i.e. bin1 with no outgoing
  // connections)

  // data for all stages
  std::vector<std::uint32_t> stageData;
  // defines the start and end of stage contents
  std::vector<std::size_t> stageOffsets;

  stageOffsets.reserve(51);  // 50 stages -> offsets needs +1
  stageOffsets.push_back(0);

  std::vector<std::uint32_t> exitBins;
  while (!binMap.empty()) {
    exitBins.clear();

    // 2a. find all bins with zero outgoing links

    for (const auto& bl : binMap) {
      auto& binLinks = bl.second;
      auto& outLinks = binLinks.first;

      if (!outLinks.empty()) {
        continue;
      }

      exitBins.push_back(bl.first);
    }

    // order exit bins (important for graph building)
    std::sort(exitBins.begin(), exitBins.end());

    // 2b. add a new stage: vector of bin1

    stageData.insert(stageData.end(), exitBins.begin(), exitBins.end());
    stageOffsets.push_back(stageData.size());

    // 2c. remove links : graph ablation

    for (const std::uint32_t& bin1Key : exitBins) {
      auto p1 = binMap.find(bin1Key);
      if (p1 == binMap.end()) {
        continue;
      }
      auto& bin1Links = p1->second;

      for (const std::uint32_t bin2Key : bin1Links.second) {
        const auto p2 = binMap.find(bin2Key);
        if (p2 == binMap.end()) {
          continue;
        }

        auto& binLinks = p2->second;
        std::vector<std::uint32_t>& links = binLinks.first;

        std::erase(links, bin1Key);
      }
    }

    // 2d. finally, remove all exit bin1s from the map

    for (auto bin1Key : exitBins) {
      binMap.erase(bin1Key);
    }
  }

  // 3. Refill binGroups with staged bin pair collections.

  m_binGroups.clear();

  // number of stages:
  const std::size_t nStages = stageOffsets.size() - 1;

  // reverse order filling
  for (std::size_t stageIndex = nStages; stageIndex-- > 0;) {
    const std::size_t begin = stageOffsets[stageIndex];
    const std::size_t end = stageOffsets[stageIndex + 1];

    for (std::size_t k = begin; k < end; ++k) {
      const std::uint32_t bin1Idx = stageData[k];

      const auto p = originalBinMap.find(bin1Idx);
      if (p == originalBinMap.end()) {
        continue;
      }

      const auto& binLists = p->second;
      const std::vector<std::uint32_t>& bin2List = binLists.second;

      // store the group
      m_binGroups.push_back({bin1Idx, bin2List});
    }
  }
}

const detail::GbtsLayer* GbtsGeometry::layerById(std::uint32_t id) const {
  if (const auto it = m_layerFromUserIdMap.find(id);
      it != m_layerFromUserIdMap.end()) {
    return &m_layers.at(it->second);
  }
  return nullptr;
}

const detail::GbtsLayer& GbtsGeometry::layerByIndex(std::uint32_t idx) const {
  return m_layers.at(idx);
}

const detail::GbtsLayer& GbtsGeometry::createLayer(
    const GbtsLayerDescription& layerDescription, std::uint32_t bin0) {
  const std::uint32_t layerIndex = m_layers.size();
  detail::GbtsLayer& ref =
      m_layers.emplace_back(layerDescription, m_etaBinWidth, bin0);
  m_layerFromUserIdMap.try_emplace(layerDescription.id, layerIndex);
  return ref;
}

}  // namespace Acts::Experimental
