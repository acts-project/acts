// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/GbtsBinning.hpp"
#include "Acts/Seeding/GbtsLayerConnection.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/detail/GbtsLayer.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <map>
#include <span>
#include <vector>

namespace Acts::Experimental {

class GbtsNodeStorage;
class GbtsTrackingFilter;
class GraphBasedTrackSeeder;

/// z0 range two eta bins must share a trajectory within, which fixes the bin
/// table. A separate cut from `GraphBasedTrackSeeder::Config::minZ0`.
struct GbtsZ0Range final {
  /// Minimum z0.
  float min = -168.0f * UnitConstants::mm;
  /// Maximum z0.
  float max = 168.0f * UnitConstants::mm;
};

/// Geometry helper built from layers and connectors.
class GbtsGeometry final {
 public:
  /// Constructor
  /// @param layerDescriptions Layer descriptions for the layers
  /// @param layerConnections Pairs of layers the seeder may connect
  /// @param etaBinWidth Width of the eta bins each layer is split into
  /// @param z0Range z0 range the bin table is built against
  /// @param logger Logging instance, only used during construction
  GbtsGeometry(std::span<const GbtsLayerDescription> layerDescriptions,
               std::span<const GbtsLayerConnection> layerConnections,
               float etaBinWidth, const GbtsZ0Range& z0Range = {},
               const Logger& logger = getDummyLogger());

  /// Get the number of layers the geometry was built from
  /// @return The layer count
  std::size_t numLayers() const { return m_layers.size(); }

  /// Get the number of eta bins across all layers
  /// @return The bin count
  std::uint32_t numBins() const { return m_nEtaBins; }

  /// Get the description a layer was built from
  /// @param idx Layer index
  /// @return Reference to the layer description
  const GbtsLayerDescription& layerDescription(std::uint32_t idx) const {
    return layerByIndex(idx).layerDescription();
  }

  /// Get the eta binning the geometry gave a layer
  /// @param idx Layer index
  /// @return Reference to the layer's binning
  const GbtsLayerBinning& layerBinning(std::uint32_t idx) const {
    return layerByIndex(idx).binning();
  }

  /// Get the eta bin pairs the graph is built over, in build order
  /// @return The bin groups
  std::span<const GbtsBinGroup> binGroups() const { return m_binGroups; }

 private:
  // The layer object is shared only with the classes that build the graph.
  friend class GbtsNodeStorage;
  friend class GbtsTrackingFilter;
  friend class GraphBasedTrackSeeder;

  /// Get layer by ID
  /// @param id Layer ID
  /// @return Pointer to layer or nullptr
  const detail::GbtsLayer* layerById(std::uint32_t id) const;

  /// Get layer by index
  /// @param idx Layer index
  /// @return Reference to layer
  const detail::GbtsLayer& layerByIndex(std::uint32_t idx) const;

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
  std::vector<GbtsBinGroup> m_binGroups;
};

}  // namespace Acts::Experimental
