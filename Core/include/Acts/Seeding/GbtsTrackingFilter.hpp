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

#include <array>
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

    /// Across the strips only. The along-strip resolution is not configured:
    /// it is derived per node from the node's own strip. See `update`.
    float sigmaXStrip = 0.03 * Acts::UnitConstants::mm;

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

    /// Initial variance of the transverse state (y, dy/dx, d2y/dx2).
    std::array<float, 3> initialVarianceX = {0.25f, 0.001f, 0.001f};
    /// Initial variance of the longitudinal state (z, dz/dr).
    std::array<float, 2> initialVarianceY = {1.5f, 0.001f};

    /// Calibrate a path's strip space points with the path's own direction,
    /// as the doublet stage and the seed parameter estimate do. A stored strip
    /// space point is where the two strips cross for a direction it does not
    /// know, which is up to about 20 mm off along the strip at low pT; fitted
    /// as stored, it biases the path's slope, and the filter then refuses true
    /// pixels further in.
    bool calibrateStrips = false;
    /// How far along a strip a calibrated crossing may land, as a multiple of
    /// the strip half-length, as `GbtsGraphBuilder::Config`.
    float maxStripLengthFraction = 1.1f;
    /// Search the pixels of a path inwards first, then pick up its strip space
    /// points outwards from the outermost pixel, with the state the pixels
    /// gave. The pixel chains then do not depend on a fit started in the
    /// strips, and the strips are calibrated with a direction the pixels set.
    /// Chains of strips alone are searched inwards as before.
    bool pixelsInwardStripsOutward = false;
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

  /// Configuration for the tracking filter.
  Config m_cfg{};

  /// GBTS geometry for layer information
  std::shared_ptr<const GbtsGeometry> m_geometry;

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

  /// Update an edge state in the inward frame with the edge's outer node: the
  /// state is carried outwards from the edge's inner node, as the strips are
  /// picked up with `pixelsInwardStripsOutward`.
  /// @param nodeView View of the node positions and layers
  /// @param pS Edge whose outer node to update with
  /// @param ts Edge state to update
  /// @return Success flag
  bool updateOutward(const detail::GbtsNodeView& nodeView,
                     const detail::GbtsEdge& pS,
                     detail::GbtsEdgeState& ts) const;

  /// Update an edge state with `newNode`, reached from `previousNode`.
  /// @param nodeView View of the node positions and layers
  /// @param newNode Node to update with
  /// @param previousNode Node the state was last updated with
  /// @param ts Edge state to update
  /// @return Success flag
  bool updateImpl(const detail::GbtsNodeView& nodeView,
                  SpacePointIndex newNode, SpacePointIndex previousNode,
                  detail::GbtsEdgeState& ts) const;

  /// Get layer type from layer index
  /// @param layerIndex Layer index
  /// @return Layer type
  GbtsLayerType getLayerType(GbtsLayerIndex layerIndex) const;
};

}  // namespace Acts::Experimental
