// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding2/GbtsLayerConnection.hpp"

#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>

namespace Acts::Experimental {

/// GBTS layer types
enum class GbtsLayerType { Barrel = 0, Endcap = 1 };

/// Lightweight layer description for GBTS geometry.
struct GbtsLayerDescription final {
  /// Combined subdetector ID.
  std::int32_t id{};
  /// Layer type (barrel or endcap).
  GbtsLayerType type{};
  /// Reference coordinate (z for barrel, r for endcap).
  float refCoord{};
  /// Minimum boundary coordinate.
  float minBound{};
  /// Maximum boundary coordinate.
  float maxBound{};
};

/// Layer helper with eta-bin access for GBTS seeding.
class GbtsLayer final {
 public:
  /// @param layerDescription Layer description for the layer
  /// @param etaBinWidth Eta bin width
  /// @param bin0 Starting bin index
  GbtsLayer(const GbtsLayerDescription& layerDescription, float etaBinWidth,
            std::int32_t bin0);

  /// Get eta bin for given z and r coordinates
  /// @param zh Z coordinate
  /// @param rh Radius coordinate
  /// @return Eta bin index
  std::int32_t getEtaBin(float zh, float rh) const;

  /// Get number of bins
  /// @return Number of bins
  std::int32_t numOfBins() const {
    return static_cast<std::int32_t>(m_bins.size());
  }

  /// Get bins
  /// @return Vector of bin indices
  const std::vector<std::int32_t>& bins() const { return m_bins; }

  /// Get the layer description
  /// @return Reference to the layer description
  const GbtsLayerDescription& layerDescription() const {
    return m_layerDescription;
  }

  /// Verify bin compatibility
  /// @param otherLayer Other layer to compare with
  /// @param b1 First bin index
  /// @param b2 Second bin index
  /// @param minZ0 Minimum z0 coordinate
  /// @param maxZ0 Maximum z0 coordinate
  /// @return True if bins are compatible
  bool checkCompatibility(const GbtsLayer& otherLayer, std::uint32_t b1,
                          std::uint32_t b2, float minZ0, float maxZ0) const;

 private:
  /// Layer description
  GbtsLayerDescription m_layerDescription;

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

/// Geometry helper built from layers and connectors.
class GbtsGeometry final {
  // map key is a bin
  // pair corresponds to outgoing and incoming bins that the current bin can
  // connect to
  using BinConnections =
      std::unordered_map<std::uint32_t, std::pair<std::vector<std::uint32_t>,
                                                  std::vector<std::uint32_t>>>;

 public:
  /// Constructor
  /// @param layerDescriptions Layer descriptions for the layers
  /// @param layerConnections Layer connections map
  GbtsGeometry(const std::vector<GbtsLayerDescription>& layerDescriptions,
               const GbtsLayerConnectionMap& layerConnections);

  /// Get number of eta bins
  /// @return Number of eta bins
  std::uint32_t numBins() const { return m_nEtaBins; }

  /// Get number of layers
  /// @return Number of layers
  std::uint32_t numLayers() const {
    return static_cast<std::uint32_t>(m_layers.size());
  }

  /// Get bin groups
  /// @return Bin groups vector
  const std::vector<std::pair<std::uint32_t, std::vector<std::uint32_t>>>&
  binGroups() const {
    return m_binGroups;
  }

  /// Get layer by ID
  /// @param id Layer ID
  /// @return Pointer to layer or nullptr
  const GbtsLayer* layerById(std::uint32_t id) const;

  /// Get layer by index
  /// @param idx Layer index
  /// @return Reference to layer
  const GbtsLayer& layerByIndex(std::int32_t idx) const;

  /// Get layer ID by index
  /// @param idx Layer index
  /// @return Layer ID
  inline std::uint32_t layerIdByIndex(std::uint32_t idx) const {
    return m_layers.at(idx).layerDescription().id;
  }

 private:
  /// @param layerDescription Layer description for the layer
  /// @param bin0 Starting bin index
  /// @return Reference to the newly added layer
  const GbtsLayer& createLayer(const GbtsLayerDescription& layerDescription,
                               std::uint32_t bin0);

  /// Eta bin width
  float m_etaBinWidth{};

  /// Layer array
  std::vector<GbtsLayer> m_layers;
  /// Layer per user ID map
  std::map<std::uint32_t, std::uint32_t> m_layerFromUserIdMap;
  /// Number of eta bins
  std::uint32_t m_nEtaBins{};

  /// Bin groups
  std::vector<std::pair<std::uint32_t, std::vector<std::uint32_t>>> m_binGroups;
};

}  // namespace Acts::Experimental
