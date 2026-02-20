// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/TrackFinding/GbtsConnector.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <vector>

namespace Acts::Experimental {

/// Lightweight silicon layer description for GBTs geometry.
struct TrigInDetSiLayer final {
  /// Constructor.
  /// @param subdet_ Subdetector identifier
  /// @param type_ Layer type (0: barrel, +/-n: endcap)
  /// @param center_ Reference coordinate
  /// @param min_ Minimum boundary
  /// @param max_ Maximum boundary
  TrigInDetSiLayer(std::int32_t subdet_, std::int16_t type_, float center_,
                   float min_, float max_)
      : subdet(subdet_),
        type(type_),
        refCoord(center_),
        minBound(min_),
        maxBound(max_) {}

  /// Combined subdetector ID.
  std::int32_t subdet{};  // combined ID
  /// Layer type (0: barrel, +/-n: endcap).
  std::int32_t type{};  // 0: barrel, +/-n : endcap
  /// Reference coordinate (z for barrel, r for endcap).
  float refCoord{};
  /// Minimum boundary coordinate.
  float minBound{};
  /// Maximum boundary coordinate.
  float maxBound{};
};

/// Layer helper with eta-bin access for GBTs seeding.
class GbtsLayer final {
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
  /// @return Reference to the silicon layer
  const TrigInDetSiLayer& getLayer() const { return *m_layer; }

  /// Verify bin compatibility
  /// @param pL Other layer to compare with
  /// @param b1 First bin index
  /// @param b2 Second bin index
  /// @param minZ0 Minimum z0 coordinate
  /// @param maxZ0 Maximum z0 coordinate
  /// @return True if bins are compatible
  bool verifyBin(const GbtsLayer& pL, std::uint32_t b1, std::uint32_t b2,
                 float minZ0, float maxZ0) const;

 private:
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
class GbtsGeometry final {
  // map key is a bin
  // pair corresponds to outgoing and incoming bins that the current bin can
  // connect to
  using BinConnections =
      std::map<std::uint32_t, std::pair<std::vector<std::uint32_t>,
                                        std::vector<std::uint32_t>>>;

 public:
  /// Constructor
  /// @param layerGeometry Silicon layers for geometry
  /// @param conn Connector for layer connections
  GbtsGeometry(const std::vector<TrigInDetSiLayer>& layerGeometry,
               const std::unique_ptr<GbtsConnector>& conn);

  /// Get layer by key
  /// @param key Layer key
  /// @return Pointer to layer or nullptr
  const GbtsLayer* getGbtsLayerByKey(std::uint32_t key) const;
  /// Get layer by index
  /// @param idx Layer index
  /// @return Reference to layer
  const GbtsLayer& getGbtsLayerByIndex(std::int32_t idx) const;

  /// Get layer key by index
  /// @param idx Layer index
  /// @return Layer key
  inline std::uint32_t getGbtsLayerKeyByIndex(std::uint32_t idx) const {
    return m_layerKeys.at(idx);
  }

  /// Get number of eta bins
  /// @return Number of eta bins
  std::uint32_t numBins() const { return m_nEtaBins; }
  /// Get number of layers
  /// @return Number of layers
  std::uint32_t numLayers() const { return m_layArray.size(); }
  /// Get bin groups
  /// @return Bin groups vector
  const std::vector<std::pair<std::uint32_t, std::vector<std::uint32_t>>>&
  binGroups() const {
    return m_binGroups;
  }

 private:
  /// Add new layer
  /// @param l Silicon layer to add
  /// @param bin0 Starting bin index
  /// @return Reference to the newly added layer
  const GbtsLayer& addNewLayer(const TrigInDetSiLayer& l, std::uint32_t bin0);

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
  std::vector<std::pair<std::uint32_t, std::vector<std::uint32_t>>> m_binGroups;
};

}  // namespace Acts::Experimental
