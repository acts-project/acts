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
#include "Acts/Seeding/detail/GbtsFilterTypes.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace Acts::Experimental {

class GbtsGeometry;
class GraphBasedTrackSeeder;

/// Tracking filter operating on the GBTS edge graph.
class GbtsTrackingFilter final {
 public:
  /// Configuration for the tracking filter.
  struct Config {
    /// Multiple scattering sigma (for 900 MeV track at eta=0).
    float sigmaMS = 0.016;
    /// Radiation length fraction per layer (2.5% per layer).
    float radLen = 0.025;

    /// Measurement uncertainty in x direction.
    float sigmaX = 0.08 * Acts::UnitConstants::mm;
    /// Measurement uncertainty in y direction.
    float sigmaY = 0.25 * Acts::UnitConstants::mm;

    /// Measurement weight in x direction.
    float weightX = 0.5;
    /// Measurement weight in y direction.
    float weightY = 0.5;

    /// Maximum delta chi2 in x direction.
    float maxDChi2X = 5.0;
    /// Maximum delta chi2 in y direction.
    float maxDChi2Y = 6.0;

    /// Chi2 penalty for adding a hit.
    float addHit = 14.0;

    /// Maximum track curvature.
    float maxCurvature = 1e-3f / Acts::UnitConstants::mm;
    /// Maximum longitudinal impact parameter.
    float maxZ0 = 170.0 * Acts::UnitConstants::mm;
  };

  /// @param config Configuration for seed finder
  /// @param geometry GBTS geometry for layer information
  /// @param logger Logging instance
  GbtsTrackingFilter(const Config& config,
                     const std::shared_ptr<const GbtsGeometry>& geometry,
                     std::unique_ptr<const Logger> logger = getDefaultLogger(
                         "GbtsTrackingFilter", Logging::Level::INFO));

 private:
  // Only the seeder walks the edge graph.
  friend class GraphBasedTrackSeeder;

  using State = detail::GbtsFilterState;

  /// Follow track starting from edge
  /// @param state Tracking filter state
  /// @param nodeView View of the node positions and layers
  /// @param sb Edge storage
  /// @param pS Starting edge
  /// @return Final edge state after following the track
  detail::GbtsEdgeState followTrack(State& state,
                                    const detail::GbtsNodeView& nodeView,
                                    std::vector<detail::GbtsEdge>& sb,
                                    detail::GbtsEdge& pS) const;

  /// Configuration for the tracking filter.
  Config m_cfg{};

  /// GBTS geometry for layer information
  std::shared_ptr<const GbtsGeometry> m_geometry;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Access the logging instance
  /// @return Logger reference
  const Logger& logger() const { return *m_logger; }

  /// Propagate edge state
  /// @param state Tracking filter state
  /// @param nodeView View of the node positions and layers
  /// @param sb Edge storage
  /// @param pS Edge to propagate from
  /// @param ts Edge state to update
  void propagate(State& state, const detail::GbtsNodeView& nodeView,
                 std::vector<detail::GbtsEdge>& sb, detail::GbtsEdge& pS,
                 detail::GbtsEdgeState& ts) const;

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
  GbtsLayerType getLayerType(std::uint32_t layerIndex) const;
};

}  // namespace Acts::Experimental
