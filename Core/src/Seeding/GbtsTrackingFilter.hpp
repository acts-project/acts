// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding/GbtsLayerDescription.hpp"
#include "Acts/Seeding/GraphBasedTrackSeeder.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

#include "GbtsFilterTypes.hpp"

namespace Acts::detail {

class GbtsGeometryImpl;
class GbtsSeederImpl;

/// Tracking filter operating on the GBTS edge graph.
class GbtsTrackingFilter final {
 public:
  /// Configuration for the tracking filter.
  using Config = Acts::GraphBasedTrackSeeder::FilterConfig;

  /// @param config Configuration for seed finder
  /// @param geometry GBTS geometry for layer information
  /// @param logger Logging instance
  GbtsTrackingFilter(const Config& config,
                     const std::shared_ptr<const GbtsGeometryImpl>& geometry,
                     std::unique_ptr<const Logger> logger = getDefaultLogger(
                         "GbtsTrackingFilter", Logging::Level::INFO));

 private:
  // Only the seeder walks the edge graph.
  friend class GbtsSeederImpl;

  /// Configuration for the tracking filter.
  Config m_cfg{};

  /// GBTS geometry for layer information
  std::shared_ptr<const GbtsGeometryImpl> m_geometry;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Access the logging instance
  /// @return Logger reference
  const Logger& logger() const { return *m_logger; }

  /// Update edge state with edge
  /// @param nodeView View of the node positions and layers
  /// @param pS Edge to update with
  /// @param ts Edge state to update
  /// @return Success flag
  bool update(const detail::GbtsNodeView& nodeView, const detail::GbtsEdge& pS,
              detail::GbtsEdgeState& ts) const;

  /// Get layer type from layer index
  /// @param layerIndex Layer index
  /// @return Layer type
  GbtsLayerType getLayerType(GbtsLayerIndex layerIndex) const;
};

}  // namespace Acts::detail
