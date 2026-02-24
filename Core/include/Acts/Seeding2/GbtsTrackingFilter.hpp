// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding2/GbtsConfig.hpp"
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

  /// Constructor
  /// @param config Configuration for seed finder
  /// @param layerGeometry Geometry layer description
  /// @param sb Edge storage
  GbtsTrackingFilter(const GbtsConfig& config,
                     const std::vector<TrigInDetSiLayer>& layerGeometry,
                     std::vector<GbtsEdge>& sb);

  /// Follow track starting from edge
  /// @param pS Starting edge
  /// @return Final edge state after following the track
  GbtsEdgeState followTrack(GbtsEdge& pS);

 private:
  /// Propagate edge state
  /// @param pS Edge to propagate from
  /// @param ts Edge state to update
  void propagate(GbtsEdge& pS, GbtsEdgeState& ts);

  /// Update edge state with edge
  /// @param pS Edge to update with
  /// @param ts Edge state to update
  /// @return Success flag
  bool update(const GbtsEdge& pS, GbtsEdgeState& ts) const;

  /// Get layer type from layer index
  /// @param l Layer index
  /// @return Layer type
  std::uint32_t getLayerType(std::uint32_t l) const;

  /// Configuration for seed finder
  const GbtsConfig* m_cfg{};

  /// Geometry layer description
  const std::vector<TrigInDetSiLayer>* m_layerGeometry{};

  /// Edge storage
  std::vector<GbtsEdge>* m_segStore{};

  /// State vector
  std::vector<GbtsEdgeState*> m_stateVec;

  /// State storage array
  std::array<GbtsEdgeState, GbtsMaxEdgeState> m_stateStore{};

  /// Global state counter
  std::uint32_t m_globalStateCounter{0};
};

}  // namespace Acts::Experimental
