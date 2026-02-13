// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsDataStorage.hpp"

#include "Acts/Seeding/GbtsGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numbers>
#include <utility>

namespace Acts::Experimental {

GbtsEtaBin::GbtsEtaBin() {
  // TODO config
  vn.reserve(1000);
}

void GbtsEtaBin::sortByPhi() {
  // TODO config
  std::array<std::vector<std::pair<float, const GbtsNode*>>, 32> phiBuckets;

  // TODO config
  const std::uint32_t nBuckets = 31;

  for (const GbtsNode* n : vn) {
    const std::uint32_t bIdx = static_cast<std::uint32_t>(
        0.5 * nBuckets *
        (n->phi / static_cast<float>(std::numbers::pi) + 1.0f));
    phiBuckets[bIdx].emplace_back(n->phi, n);
  }

  for (auto& b : phiBuckets) {
    std::ranges::sort(b);
  }

  std::uint32_t idx = 0;
  for (const auto& b : phiBuckets) {
    for (const auto& p : b) {
      vn[idx++] = p.second;
    }
  }
}

void GbtsEtaBin::initializeNodes() {
  if (vn.empty()) {
    return;
  }

  params.resize(vn.size());

  in.resize(vn.size());
  for (auto& v : in) {
    // TODO config
    v.reserve(50);  // reasonably high number of incoming edges per node
  }

  std::transform(vn.begin(), vn.end(), params.begin(), [](const GbtsNode* pN) {
    return std::array<float, 5>{-100.0, 100.0, pN->phi, pN->r, pN->z};
  });

  const auto [min_iter, max_iter] = std::ranges::minmax_element(
      vn, {}, [](const GbtsNode* s) { return s->r; });
  minRadius = (*min_iter)->r;
  maxRadius = (*max_iter)->r;
}

void GbtsEtaBin::generatePhiIndexing(float dphi) {
  for (std::uint32_t nIdx = 0; nIdx < vn.size(); nIdx++) {
    const float phi = params[nIdx][2];
    if (phi <= std::numbers::pi - dphi) {
      continue;
    }
    vPhiNodes.emplace_back(phi - 2 * std::numbers::pi, nIdx);
  }

  for (std::uint32_t nIdx = 0; nIdx < vn.size(); nIdx++) {
    const float phi = params[nIdx][2];
    vPhiNodes.emplace_back(phi, nIdx);
  }

  for (std::uint32_t nIdx = 0; nIdx < vn.size(); nIdx++) {
    const float phi = params[nIdx][2];
    if (phi >= -std::numbers::pi + dphi) {
      break;
    }
    vPhiNodes.emplace_back(phi + 2 * std::numbers::pi, nIdx);
  }
}

GbtsDataStorage::GbtsDataStorage(const GbtsConfig& config,
                                 std::shared_ptr<const GbtsGeometry> geometry,
                                 GbtsMLLookupTable mlLut)
    : m_geo(std::move(geometry)), m_cfg(config), m_mlLut(std::move(mlLut)) {
  // parse the look up table if useMl is true

  m_etaBins.resize(m_geo->numBins());
}

std::uint32_t GbtsDataStorage::loadPixelGraphNodes(
    const std::uint16_t layerIndex, const std::span<const GbtsNode> coll,
    const bool useMl) {
  std::uint32_t nLoaded = 0;

  const GbtsLayer& pL = m_geo->getGbtsLayerByIndex(layerIndex);

  const bool isBarrel = pL.getLayer().type == 0;

  for (const GbtsNode& node : coll) {
    const std::int32_t binIndex = pL.getEtaBin(node.z, node.r);

    if (binIndex == -1) {
      continue;
    }

    if (isBarrel) {
      m_etaBins.at(binIndex).vn.push_back(&node);
    } else {
      if (useMl) {
        const float clusterWidth = node.pcw;
        if (clusterWidth > m_cfg.maxEndcapClusterWidth) {
          continue;
        }
      }
      m_etaBins.at(binIndex).vn.push_back(&node);
    }

    nLoaded++;
  }

  return nLoaded;
}

std::uint32_t GbtsDataStorage::loadStripGraphNodes(
    const std::uint16_t layerIndex, const std::span<const GbtsNode> coll) {
  std::uint32_t nLoaded = 0;

  const GbtsLayer& pL = m_geo->getGbtsLayerByIndex(layerIndex);

  for (const GbtsNode& node : coll) {
    const std::int32_t binIndex = pL.getEtaBin(node.z, node.r);

    if (binIndex == -1) {
      continue;
    }

    m_etaBins.at(binIndex).vn.push_back(&node);
    nLoaded++;
  }

  return nLoaded;
}

std::uint32_t GbtsDataStorage::numberOfNodes() const {
  std::uint32_t n = 0;

  for (const auto& b : m_etaBins) {
    n += b.vn.size();
  }
  return n;
}

void GbtsDataStorage::sortByPhi() {
  for (GbtsEtaBin& b : m_etaBins) {
    b.sortByPhi();
  }
}

void GbtsDataStorage::initializeNodes(const bool useMl) {
  for (GbtsEtaBin& b : m_etaBins) {
    b.initializeNodes();
    if (!b.vn.empty()) {
      b.layerKey = m_geo->getGbtsLayerKeyByIndex((b.vn.front())->layer);
    }
  }

  if (!useMl) {
    return;
  }

  std::uint32_t nL = m_geo->numLayers();

  for (std::uint32_t layerIdx = 0; layerIdx < nL; layerIdx++) {
    const GbtsLayer& pL = m_geo->getGbtsLayerByIndex(layerIdx);

    // skip strips volumes: layers in range [1200X-1400X]
    if (pL.getLayer().subdet < 20000) {
      continue;
    }

    const bool isBarrel = pL.getLayer().type == 0;
    if (!isBarrel) {
      continue;
    }

    // adjusting cuts on |cot(theta)| using pre-trained LUT loaded from file

    const std::uint32_t lutSize = m_mlLut.size();

    const std::uint32_t nBins = pL.numOfBins();

    for (std::uint32_t b = 0; b < nBins; ++b) {  // loop over eta-bins in Layer
      GbtsEtaBin& B = m_etaBins.at(pL.getBins().at(b));

      if (B.empty()) {
        continue;
      }

      for (std::uint32_t nIdx = 0; nIdx < B.vn.size(); ++nIdx) {
        const float clusterWidth = B.vn[nIdx]->pcw;
        const float locPosY = B.vn[nIdx]->locPosY;

        const std::int32_t lutBinIdx =
            static_cast<std::uint32_t>(std::floor(20 * clusterWidth)) -
            1;  // lut bin width is 0.05 mm, check if this is
                // actually what we want with float conversion

        if (lutBinIdx >= static_cast<std::int32_t>(lutSize)) {
          continue;
        }
        if (lutBinIdx < 0) {
          continue;  // protect against negative index
        }

        const std::array<float, 5> lutBin = m_mlLut[lutBinIdx];

        const float dist2border = 10.0 - std::abs(locPosY);

        float min_tau = -100.0;
        float max_tau = 100.0;

        if (dist2border > 0.3f) {  // far enough from the edge
          min_tau = lutBin[1];
          max_tau = lutBin[2];
        } else {  // possible cluster shortening at a module edge
          min_tau = lutBin[3];
          max_tau = lutBin[4];
        }

        if (max_tau < 0) {  // insufficient training data
          max_tau = 100.0;  // use "no-cut" default
        }

        B.params[nIdx][0] = min_tau;
        B.params[nIdx][1] = max_tau;
      }
    }
  }
}

void GbtsDataStorage::generatePhiIndexing(const float dphi) {
  for (GbtsEtaBin& b : m_etaBins) {
    b.generatePhiIndexing(dphi);
  }
}

}  // namespace Acts::Experimental
