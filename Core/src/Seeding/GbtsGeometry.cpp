// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/GbtsGeometry.hpp"

#include <cmath>
#include <cstring>
#include <iostream>

namespace Acts::Experimental {

GbtsLayer::GbtsLayer(const TrigInDetSiLayer& ls, const float ew,
                     const std::int32_t bin0)
    : m_layer(&ls), m_etaBinWidth(ew) {
  if (m_layer->type == 0) {  // barrel
    m_r1 = m_layer->refCoord;
    m_r2 = m_layer->refCoord;
    m_z1 = m_layer->minBound;
    m_z2 = m_layer->maxBound;
  } else {  // endcap
    m_r1 = m_layer->minBound;
    m_r2 = m_layer->maxBound;
    m_z1 = m_layer->refCoord;
    m_z2 = m_layer->refCoord;
  }

  const float t1 = m_z1 / m_r1;
  const float eta1 = -std::log(std::sqrt(1 + t1 * t1) - t1);

  const float t2 = m_z2 / m_r2;
  const float eta2 = -std::log(std::sqrt(1 + t2 * t2) - t2);

  m_minEta = eta1;
  m_maxEta = eta2;

  if (m_maxEta < m_minEta) {
    m_minEta = eta2;
    m_maxEta = eta1;
  }

  m_maxEta += 1e-6;  // increasing them slightly to avoid range_check exceptions
  m_minEta -= 1e-6;

  const float deltaEta = m_maxEta - m_minEta;

  std::uint32_t binCounter = bin0;

  if (deltaEta < m_etaBinWidth) {
    m_nBins = 1;
    m_bins.push_back(binCounter++);
    m_etaBin = deltaEta;
    if (m_layer->type == 0) {  // barrel
      m_minRadius.push_back(m_layer->refCoord - 2.0);
      m_maxRadius.push_back(m_layer->refCoord + 2.0);
      m_minBinCoord.push_back(m_layer->minBound);
      m_maxBinCoord.push_back(m_layer->maxBound);
    } else {  // endcap
      m_minRadius.push_back(m_layer->minBound - 2.0);
      m_maxRadius.push_back(m_layer->maxBound + 2.0);
      m_minBinCoord.push_back(m_layer->minBound);
      m_maxBinCoord.push_back(m_layer->maxBound);
    }
  } else {
    const std::uint32_t nB =
        static_cast<std::uint32_t>(deltaEta / m_etaBinWidth);
    m_nBins = nB;
    if (deltaEta - m_etaBinWidth * nB > 0.5 * m_etaBinWidth) {
      m_nBins++;
    }

    m_etaBin = deltaEta / m_nBins;

    if (m_nBins == 1) {
      m_bins.push_back(binCounter++);
      if (m_layer->type == 0) {  // barrel
        m_minRadius.push_back(m_layer->refCoord - 2.0);
        m_maxRadius.push_back(m_layer->refCoord + 2.0);
        m_minBinCoord.push_back(m_layer->minBound);
        m_maxBinCoord.push_back(m_layer->maxBound);
      } else {  // endcap
        m_minRadius.push_back(m_layer->minBound - 2.0);
        m_maxRadius.push_back(m_layer->maxBound + 2.0);
        m_minBinCoord.push_back(m_layer->minBound);
        m_maxBinCoord.push_back(m_layer->maxBound);
      }
    } else {
      float eta = m_minEta + 0.5 * m_etaBin;

      for (std::uint32_t i = 1; i <= m_nBins; ++i) {
        m_bins.push_back(binCounter++);

        float e1 = eta - 0.5 * m_etaBin;
        float e2 = eta + 0.5 * m_etaBin;

        if (m_layer->type == 0) {  // barrel
          m_minRadius.push_back(m_layer->refCoord - 2.0);
          m_maxRadius.push_back(m_layer->refCoord + 2.0);
          const float z1 = m_layer->refCoord * std::sinh(e1);
          m_minBinCoord.push_back(z1);
          const float z2 = m_layer->refCoord * std::sinh(e2);
          m_maxBinCoord.push_back(z2);
        } else {  // endcap
          // for the positive endcap larger eta corresponds to smaller radius
          if (m_layer->refCoord > 0) {
            std::swap(e1, e2);
          }
          float r = m_layer->refCoord / std::sinh(e1);
          m_minBinCoord.push_back(r);
          m_minRadius.push_back(r - 2.0);
          r = m_layer->refCoord / std::sinh(e2);
          m_maxBinCoord.push_back(r);
          m_maxRadius.push_back(r + 2.0);
        }

        eta += m_etaBin;
      }
    }
  }
}

bool GbtsLayer::verifyBin(const GbtsLayer& pL, const std::uint32_t b1,
                          const std::uint32_t b2, const float minZ0,
                          const float maxZ0) const {
  const float z1min = m_minBinCoord.at(b1);
  const float z1max = m_maxBinCoord.at(b1);
  const float r1 = m_layer->refCoord;

  if (m_layer->type == 0 && pL.m_layer->type == 0) {  // barrel <- barrel
    const float tol = 5.0;

    const float min_b2 = pL.m_minBinCoord.at(b2);
    const float max_b2 = pL.m_maxBinCoord.at(b2);

    const float r2 = pL.m_layer->refCoord;

    const float A = r2 / (r2 - r1);
    const float B = r1 / (r2 - r1);

    const float z0_min = z1min * A - max_b2 * B;
    const float z0_max = z1max * A - min_b2 * B;

    if (z0_max < minZ0 - tol || z0_min > maxZ0 + tol) {
      return false;
    }

    return true;
  }

  if (m_layer->type == 0 && pL.m_layer->type != 0) {  // barrel <- endcap
    const float z2 = pL.m_layer->refCoord;
    const float r2max = pL.m_maxBinCoord.at(b2);
    float r2min = pL.m_minBinCoord.at(b2);

    if (r2max <= r1) {
      return false;
    }

    if (r2min <= r1) {
      r2min = r1 + 1e-3;
    }

    float z0_max = 0.0;
    float z0_min = 0.0;

    if (z2 > 0) {
      z0_max = (z1max * r2max - z2 * r1) / (r2max - r1);
      z0_min = (z1min * r2min - z2 * r1) / (r2min - r1);
    } else {
      z0_max = (z1max * r2min - z2 * r1) / (r2min - r1);
      z0_min = (z1min * r2max - z2 * r1) / (r2max - r1);
    }

    if (z0_max < minZ0 || z0_min > maxZ0) {
      return false;
    }
    return true;
  }

  return true;
}

std::int32_t GbtsLayer::getEtaBin(const float zh, const float rh) const {
  if (m_bins.size() == 1) {
    return m_bins.at(0);
  }

  const float t1 = zh / rh;
  const float eta = -std::log(std::sqrt(1 + t1 * t1) - t1);

  std::int32_t idx = static_cast<std::int32_t>((eta - m_minEta) / m_etaBin);
  if (idx < 0) {
    idx = 0;
  } else if (idx >= static_cast<std::int32_t>(m_bins.size())) {
    idx = static_cast<std::int32_t>(m_bins.size()) - 1;
  }

  return m_bins.at(idx);  // index in the global storage
}

GbtsGeometry::GbtsGeometry(const std::vector<TrigInDetSiLayer>& layerGeometry,
                           const std::unique_ptr<GbtsConnector>& conn) {
  // TODO configurable z0 range
  const float minZ0 = -168.0;
  const float maxZ0 = 168.0;

  m_etaBinWidth = conn->m_etaBin;

  for (const TrigInDetSiLayer& layer : layerGeometry) {
    const GbtsLayer& pL = addNewLayer(layer, m_nEtaBins);
    m_nEtaBins += pL.numOfBins();
  }

  // calculating bin tables in the connector...
  // calculate bin pairs for graph edge building

  std::int32_t lastBin1 = -1;

  for (const auto& [layer, vConn] : conn->m_connMap) {
    for (const auto& connection : vConn) {
      const std::uint32_t src = connection->m_src;  // n2 : the new connectors
      const std::uint32_t dst = connection->m_dst;  // n1

      const GbtsLayer* pL1 = getGbtsLayerByKey(dst);
      const GbtsLayer* pL2 = getGbtsLayerByKey(src);
      if (pL1 == nullptr) {
        std::cout << " skipping invalid dst layer " << dst << std::endl;
        continue;
      }
      if (pL2 == nullptr) {
        std::cout << " skipping invalid src layer " << src << std::endl;
        continue;
      }

      const std::uint32_t nSrcBins = pL2->numOfBins();
      const std::uint32_t nDstBins = pL1->numOfBins();

      connection->m_binTable.resize(nSrcBins * nDstBins, 0);
      // loop over bins in Layer 1
      for (std::uint32_t b1 = 0; b1 < nDstBins; ++b1) {
        // loop over bins in Layer 2
        for (std::uint32_t b2 = 0; b2 < nSrcBins; ++b2) {
          if (!pL1->verifyBin(*pL2, b1, b2, minZ0, maxZ0)) {
            continue;
          }
          const std::uint32_t address = b1 + b2 * nDstBins;
          connection->m_binTable.at(address) = 1;

          const std::int32_t bin1Idx = pL1->getBins().at(b1);
          const std::int32_t bin2Idx = pL2->getBins().at(b2);

          if (bin1Idx != lastBin1) {  // adding a new group
            m_binGroups.emplace_back(bin1Idx,
                                     std::vector<std::int32_t>(1, bin2Idx));
            lastBin1 = bin1Idx;
          } else {  // extend the last group
            m_binGroups.back().second.push_back(bin2Idx);
          }
        }
      }
    }
  }
}

const GbtsLayer* GbtsGeometry::getGbtsLayerByKey(std::uint32_t key) const {
  if (const auto it = m_layMap.find(key); it != m_layMap.end()) {
    return it->second;
  }
  return nullptr;
}

const GbtsLayer& GbtsGeometry::getGbtsLayerByIndex(std::int32_t idx) const {
  return *m_layArray.at(idx);
}

const GbtsLayer& GbtsGeometry::addNewLayer(const TrigInDetSiLayer& l,
                                           std::uint32_t bin0) {
  const std::uint32_t layerKey = l.subdet;
  const float ew = m_etaBinWidth;

  auto pHL = std::make_unique<GbtsLayer>(l, ew, bin0);
  GbtsLayer& ref = *pHL;

  m_layMap.try_emplace(layerKey, &ref);
  m_layArray.push_back(std::move(pHL));
  m_layerKeys.push_back(layerKey);
  return ref;
}

}  // namespace Acts::Experimental
