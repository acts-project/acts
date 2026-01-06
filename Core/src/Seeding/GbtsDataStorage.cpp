// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsDataStorage.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numbers>
namespace Acts::Experimental {

GbtsEtaBin::GbtsEtaBin() {
  m_in.clear();
  m_vn.clear();
  m_params.clear();
  m_vn.reserve(1000);
}

GbtsEtaBin::~GbtsEtaBin() {
  m_in.clear();
  m_vn.clear();
  m_params.clear();
}

void GbtsEtaBin::sortByPhi() {
  std::vector<std::pair<float, const GbtsNode*>> phiBuckets[32];

  int nBuckets = 31;

  for (const auto& n : m_vn) {
    int bIdx = static_cast<int>(
        0.5 * nBuckets *
        (n->phi() / static_cast<float>(std::numbers::pi) + 1.0f));
    phiBuckets[bIdx].push_back(std::make_pair(n->phi(), n));
  }

  for (auto& b : phiBuckets) {
    std::sort(b.begin(), b.end());
  }

  int idx = 0;
  for (const auto& b : phiBuckets) {
    for (const auto& p : b) {
      m_vn[idx++] = p.second;
    }
  }
}

void GbtsEtaBin::initializeNodes() {
  if (m_vn.empty()) {
    return;
  }

  m_params.resize(m_vn.size());

  m_in.resize(m_vn.size());
  for (auto& v : m_in) {
    v.reserve(50);  // reasonably high number of incoming edges per node
  }

  std::transform(
      m_vn.begin(), m_vn.end(), m_params.begin(), [](const GbtsNode* pN) {
        std::array<float, 5> a = {-100.0, 100.0, pN->phi(), pN->r(), pN->z()};
        return a;
      });

  auto [min_iter, max_iter] = std::minmax_element(
      m_vn.begin(), m_vn.end(),
      [](const GbtsNode* s, const GbtsNode* s1) { return (s->r() < s1->r()); });
  m_maxRadius = (*max_iter)->r();
  m_minRadius = (*min_iter)->r();
}

void GbtsEtaBin::generatePhiIndexing(float dphi) {
  for (unsigned int nIdx = 0; nIdx < m_vn.size(); nIdx++) {
    float phi = m_params[nIdx][2];
    if (phi <= std::numbers::pi - dphi) {
      continue;
    }
    m_vPhiNodes.push_back(
        std::pair<float, unsigned int>(phi - 2 * std::numbers::pi, nIdx));
  }

  for (unsigned int nIdx = 0; nIdx < m_vn.size(); nIdx++) {
    float phi = m_params[nIdx][2];
    m_vPhiNodes.push_back(std::pair<float, unsigned int>(phi, nIdx));
  }

  for (unsigned int nIdx = 0; nIdx < m_vn.size(); nIdx++) {
    float phi = m_params[nIdx][2];
    if (phi >= -std::numbers::pi + dphi) {
      break;
    }
    m_vPhiNodes.push_back(
        std::pair<float, unsigned int>(phi + 2 * std::numbers::pi, nIdx));
  }
}
GbtsDataStorage::GbtsDataStorage(
    const GbtsGeometry& geometry, const SeedFinderGbtsConfig& config,
    const std::vector<std::array<float, 5>>& parsedLutFile)
    : m_geo(geometry), m_config(config), m_mlLUT(parsedLutFile) {
  // parse the look up table if useML is true

  m_etaBins.resize(geometry.num_bins());
}

GbtsDataStorage::~GbtsDataStorage() = default;

int GbtsDataStorage::loadPixelGraphNodes(short layerIndex,
                                         const std::vector<GbtsNode>& coll,
                                         bool useML) {
  int nLoaded = 0;

  const GbtsLayer* pL = m_geo.getGbtsLayerByIndex(layerIndex);

  if (pL == nullptr) {
    return -1;
  }

  bool isBarrel = (pL->m_layer.m_type == 0);

  for (const auto& node : coll) {
    int binIndex = pL->getEtaBin(node.z(), node.r());

    if (binIndex == -1) {
      continue;
    }

    if (isBarrel) {
      m_etaBins.at(binIndex).m_vn.push_back(&node);
    } else {
      if (useML) {
        float cluster_width = node.pixelClusterWidth();
        if (cluster_width > m_config.max_endcap_clusterwidth) {
          continue;
        }
      }
      m_etaBins.at(binIndex).m_vn.push_back(&node);
    }

    nLoaded++;
  }

  return nLoaded;
}

int GbtsDataStorage::loadStripGraphNodes(short layerIndex,
                                         const std::vector<GbtsNode>& coll) {
  int nLoaded = 0;

  const GbtsLayer* pL = m_geo.getGbtsLayerByIndex(layerIndex);

  if (pL == nullptr) {
    return -1;
  }

  for (const auto& node : coll) {
    int binIndex = pL->getEtaBin(node.z(), node.r());

    if (binIndex == -1) {
      continue;
    }

    m_etaBins.at(binIndex).m_vn.push_back(&node);
    nLoaded++;
  }

  return nLoaded;
}

unsigned int GbtsDataStorage::numberOfNodes() const {
  unsigned int n = 0;

  for (const auto& b : m_etaBins) {
    n += b.m_vn.size();
  }
  return n;
}

void GbtsDataStorage::sortByPhi() {
  for (auto& b : m_etaBins) {
    b.sortByPhi();
  }
}

void GbtsDataStorage::initializeNodes(bool useML) {
  for (auto& b : m_etaBins) {
    b.initializeNodes();
    if (!b.m_vn.empty()) {
      b.m_layerKey = m_geo.getGbtsLayerKeyByIndex((*b.m_vn.begin())->m_layer);
    }
  }

  if (!useML) {
    return;
  }

  unsigned int nL = m_geo.num_layers();

  for (unsigned int layerIdx = 0; layerIdx < nL; layerIdx++) {
    const GbtsLayer* pL = m_geo.getGbtsLayerByIndex(layerIdx);

    if (pL->m_layer.m_subdet <
        20000) {  // skip strips volumes: layers in range [1200X-1400X]
      continue;
    }

    bool isBarrel = (pL->m_layer.m_type == 0);

    if (!isBarrel) {
      continue;
    }

    // adjusting cuts on |cot(theta)| using pre-trained LUT loaded from file

    int lutSize = m_mlLUT.size();

    int nBins = pL->m_bins.size();

    for (int b = 0; b < nBins; b++) {  // loop over eta-bins in Layer

      GbtsEtaBin& B = m_etaBins.at(pL->m_bins.at(b));

      if (B.empty()) {
        continue;
      }

      for (unsigned int nIdx = 0; nIdx < B.m_vn.size(); nIdx++) {
        float cluster_width = B.m_vn[nIdx]->pixelClusterWidth();
        float locPosY = B.m_vn[nIdx]->localPositionY();

        int lutBinIdx = static_cast<int>(std::floor(20 * cluster_width)) -
                        1;  // lut bin width is 0.05 mm, check if this is
                            // actually what we want with float conversion

        if (lutBinIdx >= lutSize) {
          continue;
        }
        if (lutBinIdx < 0) {
          continue;  // protect against negative index
        }

        const std::array<float, 5> lutBin = m_mlLUT[lutBinIdx];

        float dist2border = 10.0 - std::abs(locPosY);

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

        B.m_params[nIdx][0] = min_tau;
        B.m_params[nIdx][1] = max_tau;
      }
    }
  }
}

void GbtsDataStorage::generatePhiIndexing(float dphi) {
  for (auto& b : m_etaBins) {
    b.generatePhiIndexing(dphi);
  }
}

}  // namespace Acts::Experimental
