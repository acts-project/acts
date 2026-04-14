// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Seeding2/GbtsDataStorage.hpp"
#include "Acts/Seeding2/GbtsGeometry.hpp"

#include <cstring>
#include <vector>

namespace Acts::Experimental {

/// Per-edge tracking state used by the GBTS filter.
struct GbtsEdgeState final {
 public:
  GbtsEdgeState() = default;

  /// Constructor with initialization flag
  /// @param f Initialization flag
  explicit GbtsEdgeState(bool f) : initialized(f) {}

  /// Initialize from edge
  /// @param pS Edge to initialize from
  void initialize(const GbtsEdge& pS);

  /// Initialization flag
  bool initialized{false};

  /// Score for comparison
  float j{};

  /// Vector of edges in the track
  std::vector<GbtsEdge*> vs;

  /// State vector X
  std::array<float, 3> x{};
  /// State vector Y
  std::array<float, 2> y{};
  /// Covariance matrix for X
  std::array<std::array<float, 3>, 3> cx{};
  /// Covariance matrix for Y
  std::array<std::array<float, 2>, 2> cy{};
  /// Reference x coordinate
  float refX{};
  /// Reference y coordinate
  float refY{};
  /// Cosine of rotation angle
  float c{};
  /// Sine of rotation angle
  float s{};
};

/// Tracking filter operating on the GBTS edge graph.
class GbtsTrackingFilter final {
 public:
  /// Maximum number of edge states
  static constexpr std::uint32_t GbtsMaxEdgeState = 2500;

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

  /// State for the tracking filter, containing edge states and a global
  /// counter.
  struct State {
    /// State vector
    std::vector<GbtsEdgeState*> stateVec;

    /// State storage array
    std::array<GbtsEdgeState, GbtsMaxEdgeState> stateStore{};

    /// Global state counter
    std::uint32_t globalStateCounter{0};
  };

  /// @param config Configuration for seed finder
  /// @param geometry GBTS geometry for layer information
  GbtsTrackingFilter(const Config& config,
                     const std::shared_ptr<const GbtsGeometry>& geometry);

  /// Follow track starting from edge
  /// @param state Tracking filter state
  /// @param sb Edge storage
  /// @param pS Starting edge
  /// @return Final edge state after following the track
  GbtsEdgeState followTrack(State& state, std::vector<GbtsEdge>& sb,
                            GbtsEdge& pS) const;

 private:
  /// Configuration for the tracking filter.
  Config m_cfg{};

  /// GBTS geometry for layer information
  std::shared_ptr<const GbtsGeometry> m_geometry;

  /// Propagate edge state
  /// @param state Tracking filter state
  /// @param sb Edge storage
  /// @param pS Edge to propagate from
  /// @param ts Edge state to update
  void propagate(State& state, std::vector<GbtsEdge>& sb, GbtsEdge& pS,
                 GbtsEdgeState& ts) const;

  /// Update edge state with edge
  /// @param pS Edge to update with
  /// @param ts Edge state to update
  /// @return Success flag
  bool update(const GbtsEdge& pS, GbtsEdgeState& ts) const;

  /// Get layer type from layer index
  /// @param layerIndex Layer index
  /// @return Layer type
  GbtsLayerType getLayerType(std::uint32_t layerIndex) const;
};

}  // namespace Acts::Experimental
