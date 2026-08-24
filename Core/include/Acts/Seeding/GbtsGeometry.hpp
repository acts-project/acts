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
