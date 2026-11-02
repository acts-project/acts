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
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>

namespace Acts {

namespace detail {
class GbtsGeometryImpl;
}
class GraphBasedTrackSeeder;

/// z0 range two eta bins must share a trajectory within, which fixes the bin
/// table. A separate cut from `GraphBasedTrackSeeder::GraphConfig::minZ0`.
struct GbtsZ0Range final {
  /// Minimum z0.
  float min = -168.0f * UnitConstants::mm;
  /// Maximum z0.
  float max = 168.0f * UnitConstants::mm;
  /// Margin added to both ends of the range, for hits off the layer line.
  float z0Tolerance = 5.0f * UnitConstants::mm;
};

/// Geometry helper built from layers and connectors.
/// Copies share immutable geometry data.
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
  std::size_t numLayers() const;

  /// Get the number of eta bins across all layers
  /// @return The bin count
  std::uint32_t numBins() const;

  /// Resolve a layer id into the index GbtsNodeStorage::insert takes.
  /// @param id Layer id, as the layer descriptions carry it
  /// @return The layer's index, or nullopt if this geometry has no such layer
  std::optional<GbtsLayerIndex> layerIndex(GbtsExperimentLayerId id) const;

  /// Get the description a layer was built from
  /// @param idx Layer index
  /// @return Reference to the layer description
  const GbtsLayerDescription& layerDescription(GbtsLayerIndex idx) const;

  /// Get the eta binning the geometry gave a layer
  /// @param idx Layer index
  /// @return Reference to the layer's binning
  const GbtsLayerBinning& layerBinning(GbtsLayerIndex idx) const;

  /// Get the eta bin pairs the graph is built over, in build order
  /// @return The bin groups
  std::span<const GbtsBinGroup> binGroups() const;

 private:
  friend class GraphBasedTrackSeeder;
  std::shared_ptr<const detail::GbtsGeometryImpl> m_impl;
};

}  // namespace Acts
