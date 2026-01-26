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
  int m_subdet;  // combined ID
  int m_type;    // 0: barrel, +/-n : endcap
  float m_refCoord;
  float m_minBound, m_maxBound;

  TrigInDetSiLayer(int subdet, short int type, float center, float min,
                   float max)
      : m_subdet(subdet),
        m_type(type),
        m_refCoord(center),
        m_minBound(min),
        m_maxBound(max) {}
};

class GbtsLayer {
 public:
  GbtsLayer(const TrigInDetSiLayer& ls, float ew, std::int32_t bin0);

  // accessors
  int getEtaBin(float zh, float rh) const;

  float getMinBinRadius(int idx) const;

  float getMaxBinRadius(int idx) const;

  std::int32_t numOfBins() const { return m_bins.size(); }

  const std::vector<std::int32_t>& getBins() const { return m_bins; }

  const TrigInDetSiLayer* getLayer() const { return m_layer; }

  bool verifyBin(const GbtsLayer* pL, int b1, int b2, float min_z0,
                 float max_z0) const;

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
  std::int32_t m_nBins{};
};

class GbtsGeometry {
 public:
  GbtsGeometry(const std::vector<TrigInDetSiLayer>& layers,
               const std::unique_ptr<GbtsConnector>& conn);

  const GbtsLayer* getGbtsLayerByKey(std::uint32_t key) const;
  const GbtsLayer* getGbtsLayerByIndex(std::int32_t idx) const;

  inline unsigned int getGbtsLayerKeyByIndex(int idx) const {
    return m_layerKeys[idx];
  }

  int num_bins() const { return m_nEtaBins; }
  unsigned int num_layers() const { return m_layArray.size(); }
  const std::vector<std::pair<int, std::vector<int>>>& bin_groups() const {
    return m_binGroups;
  }

 protected:
  const GbtsLayer* addNewLayer(const TrigInDetSiLayer& l, int bin0);

  float m_etaBinWidth{};

  std::map<unsigned int, GbtsLayer*> m_layMap;
  std::vector<std::unique_ptr<GbtsLayer>> m_layArray;
  std::vector<unsigned int> m_layerKeys;
  int m_nEtaBins{};

  std::vector<std::pair<int, std::vector<int>>> m_binGroups;
};

}  // namespace Acts::Experimental
