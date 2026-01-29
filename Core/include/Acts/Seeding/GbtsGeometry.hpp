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
/// Lightweight silicon layer description for GBTs geometry.
class TrigInDetSiLayer {
 public:
  /// Combined subdetector ID.
  std::int32_t m_subdet;  // combined ID
  /// Layer type (0: barrel, +/-n: endcap).
  std::int32_t m_type;  // 0: barrel, +/-n : endcap
  /// Reference coordinate (z for barrel, r for endcap).
  float m_refCoord;
  /// Minimum boundary coordinate.
  float m_minBound;
  /// Maximum boundary coordinate.
  float m_maxBound;

  /// Constructor.
  /// @param subdet Subdetector identifier
  /// @param type Layer type (0: barrel, +/-n: endcap)
  /// @param center Reference coordinate
  /// @param min Minimum boundary
  /// @param max Maximum boundary
  TrigInDetSiLayer(std::int32_t subdet, std::int16_t type, float center,
                   float min, float max)
      : m_subdet(subdet),
        m_type(type),
        m_refCoord(center),
        m_minBound(min),
        m_maxBound(max) {}
};

/// Layer helper with eta-bin access for GBTs seeding.
class GbtsLayer {
 public:
  /// Constructor
  /// @param ls Silicon layer descriptor
  /// @param ew Eta bin width
  /// @param bin0 Starting bin index
  GbtsLayer(const TrigInDetSiLayer& ls, float ew, std::int32_t bin0);

  // accessors
  /// Get eta bin for given z and r coordinates
  /// @param zh Z coordinate
  /// @param rh Radius coordinate
  /// @return Eta bin index
  std::int32_t getEtaBin(float zh, float rh) const;

  /// Get minimum radius for bin index
  /// @param idx Bin index
  /// @return Minimum radius for the bin
  float getMinBinRadius(std::int32_t idx) const;
  /// Get maximum radius for bin index
  /// @param idx Bin index
  /// @return Maximum radius for the bin
  float getMaxBinRadius(std::int32_t idx) const;

  /// Get number of bins
  /// @return Number of bins
  std::int32_t numOfBins() const { return m_bins.size(); }

  /// Get bins
  /// @return Vector of bin indices
  const std::vector<std::int32_t>& getBins() const { return m_bins; }

  /// Get layer
  /// @return Pointer to the silicon layer
  const TrigInDetSiLayer* getLayer() const { return m_layer; }

  /// Verify bin compatibility
  /// @param pL Other layer to compare with
  /// @param b1 First bin index
  /// @param b2 Second bin index
  /// @param min_z0 Minimum z0 coordinate
  /// @param max_z0 Maximum z0 coordinate
  /// @return True if bins are compatible
  bool verifyBin(const GbtsLayer* pL, std::uint32_t b1, std::uint32_t b2,
                 float minZ0, float maxZ0) const;

 protected:
  /// Layer
  const TrigInDetSiLayer* m_layer{};
  /// Eta-bin indices
  std::vector<std::int32_t> m_bins;
  /// Minimum radius per bin
  std::vector<float> m_minRadius;
  /// Maximum radius per bin
  std::vector<float> m_maxRadius;
  /// Minimum bin coordinate
  std::vector<float> m_minBinCoord;
  /// Maximum bin coordinate
  std::vector<float> m_maxBinCoord;

  /// Minimum eta
  float m_minEta{};
  /// Maximum eta
  float m_maxEta{};
  /// Eta bin
  float m_etaBin{};
  /// Eta bin width
  float m_etaBinWidth{};
  /// First radius coordinate
  float m_r1{};
  /// First z coordinate
  float m_z1{};
  /// Second radius coordinate
  float m_r2{};
  /// Second z coordinate
  float m_z2{};
  /// Number of bins
  std::uint32_t m_nBins{};
};

/// Geometry helper built from silicon layers and connectors.
class GbtsGeometry {
 public:
  /// Constructor
  /// @param layers Silicon layers for geometry
  /// @param conn Connector for layer connections
  GbtsGeometry(const std::vector<TrigInDetSiLayer>& layers,
               const std::unique_ptr<GbtsConnector>& conn);

  /// Get layer by key
  /// @param key Layer key
  /// @return Pointer to layer or nullptr
  const GbtsLayer* getGbtsLayerByKey(std::uint32_t key) const;
  /// Get layer by index
  /// @param idx Layer index
  /// @return Pointer to layer or nullptr
  const GbtsLayer* getGbtsLayerByIndex(std::int32_t idx) const;

  /// Get layer key by index
  /// @param idx Layer index
  /// @return Layer key
  inline std::uint32_t getGbtsLayerKeyByIndex(std::uint32_t idx) const {
    return m_layerKeys[idx];
  }

  /// Get number of eta bins
  /// @return Number of eta bins
  std::uint32_t numBins() const { return m_nEtaBins; }
  /// Get number of layers
  /// @return Number of layers
  std::uint32_t numLayers() const { return m_layArray.size(); }
  /// Get bin groups
  /// @return Bin groups vector
  const std::vector<std::pair<std::int32_t, std::vector<std::int32_t>>>&
  binGroups() const {
    return m_binGroups;
  }

 protected:
  /// Add new layer
  /// @param l Silicon layer to add
  /// @param bin0 Starting bin index
  /// @return Pointer to the newly added layer
  const GbtsLayer* addNewLayer(const TrigInDetSiLayer& l, std::uint32_t bin0);

  /// Eta bin width
  float m_etaBinWidth{};

  /// Layer map
  std::map<std::uint32_t, GbtsLayer*> m_layMap;
  /// Layer array
  std::vector<std::unique_ptr<GbtsLayer>> m_layArray;
  /// Layer keys
  std::vector<std::uint32_t> m_layerKeys;
  /// Number of eta bins
  std::uint32_t m_nEtaBins{};

  /// Bin groups
  std::vector<std::pair<std::int32_t, std::vector<std::int32_t>>> m_binGroups;
};

}  // namespace Acts::Experimental
