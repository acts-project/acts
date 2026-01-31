// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
#include "Acts/TrackFinding/GbtsConnector.hpp"

#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <vector>

namespace Acts::Experimental {
class TrigInDetSiLayer {
 public:
  std::int32_t m_subdet;  // combined ID
  std::int32_t m_type;    // 0: barrel, +/-n : endcap
  float m_refCoord;
  float m_minBound, m_maxBound;

  TrigInDetSiLayer(std::int32_t subdet, std::int16_t type, float center,
                   float min, float max)
      : m_subdet(subdet),
        m_type(type),
        m_refCoord(center),
        m_minBound(min),
        m_maxBound(max) {}
};

class GbtsLayer {
 public:
  GbtsLayer(const TrigInDetSiLayer& ls, float ew, std::int32_t bin0);

  std::int32_t getEtaBin(float zh, float rh) const;

  std::int32_t numOfBins() const { return m_bins.size(); }

  const std::vector<std::int32_t>& getBins() const { return m_bins; }

  const TrigInDetSiLayer* getLayer() const { return m_layer; }

  bool verifyBin(const GbtsLayer* pL, std::uint32_t b1, std::uint32_t b2,
                 float minZ0, float maxZ0) const;

 protected:
  const TrigInDetSiLayer* m_layer{};
  std::vector<std::int32_t> m_bins;  // eta-bin indices
  std::vector<float> m_minRadius;
  std::vector<float> m_maxRadius;
  std::vector<float> m_minBinCoord;
  std::vector<float> m_maxBinCoord;

  float m_minEta{};
  float m_maxEta{};
  float m_etaBin{};
  float m_etaBinWidth{};
  float m_r1{};
  float m_z1{};
  float m_r2{};
  float m_z2{};
  std::uint32_t m_nBins{};
};

class GbtsGeometry {
 public:
  GbtsGeometry(const std::vector<TrigInDetSiLayer>& layers,
               const std::unique_ptr<GbtsConnector>& conn);

  const GbtsLayer* getGbtsLayerByKey(std::uint32_t key) const;
  const GbtsLayer* getGbtsLayerByIndex(std::int32_t idx) const;

  inline std::uint32_t getGbtsLayerKeyByIndex(std::uint32_t idx) const {
    return m_layerKeys[idx];
  }

  std::uint32_t numBins() const { return m_nEtaBins; }
  std::uint32_t numLayers() const { return m_layArray.size(); }
  const std::vector<std::pair<std::int32_t, std::vector<std::int32_t>>>&
  binGroups() const {
    return m_binGroups;
  }

 protected:
  const GbtsLayer* addNewLayer(const TrigInDetSiLayer& l, std::uint32_t bin0);

  float m_etaBinWidth{};

  std::map<std::uint32_t, GbtsLayer*> m_layMap;
  std::vector<std::unique_ptr<GbtsLayer>> m_layArray;
  std::vector<std::uint32_t> m_layerKeys;
  std::uint32_t m_nEtaBins{};

  std::vector<std::pair<std::int32_t, std::vector<std::int32_t>>> m_binGroups;
};

}  // namespace Acts::Experimental
