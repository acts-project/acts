// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/detail/GbtsLayer.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <map>
#include <vector>

namespace Acts::Experimental {

class GbtsNodeStorage;
class GbtsTrackingFilter;
class GraphBasedTrackSeeder;

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
  /// Half width in reference direction
  float halfRefWidth = {};
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
  std::int32_t numOfBins() const { return m_bins.size(); }

  /// Get bins
  /// @return Vector of bin indices
  const std::vector<std::int32_t>& bins() const { return m_bins; }

  /// Get the layer description
  /// @return Reference to the layer description
  const GbtsLayerDescription& layerDescription() const {
    return m_layerDescription;
  }

  const float minEta() const { return m_minEta; }
  const float etaBin() const { return m_etaBin; }

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
 public:
  /// Constructor
  /// @param layerDescriptions Layer descriptions for the layers
  /// @param layerConnections Layer connections map
  /// @param logger Logging instance, only used during construction
  GbtsGeometry(const std::vector<GbtsLayerDescription>& layerDescriptions,
               const GbtsLayerConnectionMap& layerConnections,
               const Logger& logger = getDummyLogger());

 private:
  // The layer binning is shared only with the classes that build the graph.
  friend class GbtsNodeStorage;
  friend class GbtsTrackingFilter;
  friend class GraphBasedTrackSeeder;

  /// Get number of eta bins
  /// @return Number of eta bins
  std::uint32_t numBins() const { return m_nEtaBins; }

  /// Get bin groups
  /// @return Bin groups vector
  const std::vector<std::pair<std::uint32_t, std::vector<std::uint32_t>>>&
  binGroups() const {
    return m_binGroups;
  }

  /// Get layer by ID
  /// @param id Layer ID
  /// @return Pointer to layer or nullptr
  const detail::GbtsLayer* layerById(std::uint32_t id) const;

  /// Get layer by index
  /// @param idx Layer index
  /// @return Reference to layer
  const detail::GbtsLayer& layerByIndex(std::int32_t idx) const;

  /// Get layer ID by index
  /// @param idx Layer index
  /// @return Layer ID
  std::uint32_t layerIdByIndex(std::uint32_t idx) const {
    return m_layers.at(idx).layerDescription().id;
  }

  /// @param layerDescription Layer description for the layer
  /// @param bin0 Starting bin index
  /// @return Reference to the newly added layer
  const detail::GbtsLayer& createLayer(
      const GbtsLayerDescription& layerDescription, std::uint32_t bin0);

  /// Eta bin width
  float m_etaBinWidth{};

  /// Layer array
  std::vector<detail::GbtsLayer> m_layers;
  /// Layer per user ID map
  std::map<std::uint32_t, std::uint32_t> m_layerFromUserIdMap;
  /// Number of eta bins
  std::uint32_t m_nEtaBins{};

  /// Bin groups
  std::vector<std::pair<std::uint32_t, std::vector<std::uint32_t>>> m_binGroups;
};

}  // namespace Acts::Experimental
