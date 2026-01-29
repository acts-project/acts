// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
#include "Acts/Seeding/GbtsDataStorage.hpp"  //includes geo which has trigindetsilayer, may move this to trigbase
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"

#include <cmath>
#include <cstring>
#include <vector>

namespace Acts::Experimental {

/// Per-edge tracking state used by the GBTs filter.
struct GbtsEdgeState {
 public:
  /// Comparator sorting edge states by descending score.
  struct Compare {
    /// Compare two edge states by score
    /// @param s1 First edge state
    /// @param s2 Second edge state
    /// @return True if first state has higher score
    bool operator()(const struct GbtsEdgeState* s1,
                    const struct GbtsEdgeState* s2) {
      return s1->m_J > s2->m_J;
    }
  };

  GbtsEdgeState() = default;

  /// Constructor with initialization flag
  /// @param f Initialization flag
  explicit GbtsEdgeState(bool f) : m_initialized(f) {}

  /// Initialize from edge
  /// @param pS Edge to initialize from
  void initialize(GbtsEdge& pS);
  /// Clone from another state
  /// @param st Source state to clone from
  void clone(const struct GbtsEdgeState& st);

  /// Score for comparison
  /// Score for comparison
  float m_J{};

  /// Vector of edges in the track
  std::vector<GbtsEdge*> m_vs;

  /// State vector X
  float m_X[3]{};
  /// State vector Y
  float m_Y[2]{};
  /// Covariance matrix for X
  float m_Cx[3][3]{};
  /// Covariance matrix for Y
  float m_Cy[2][2]{};
  /// Reference x coordinate
  float m_refX{};
  /// Reference y coordinate
  float m_refY{};
  /// Cosine of rotation angle
  float m_c{};
  /// Sine of rotation angle
  float m_s{};

  /// Initialization flag
  bool m_initialized{false};
};

/// Tracking filter operating on the GBTs edge graph.
class GbtsTrackingFilter {
 public:
  /// Maximum number of edge states
  static constexpr std::uint32_t GbtsMaxEdgeState = 2500;
  /// Constructor
  /// @param g Geometry layer description
  /// @param sb Edge storage
  /// @param config Configuration for seed finder
  GbtsTrackingFilter(const std::vector<TrigInDetSiLayer>& g,
                     std::vector<GbtsEdge>& sb,
                     const SeedFinderGbtsConfig& config);

  /// Follow track starting from edge
  /// @param pS Starting edge
  /// @param output Output edge state
  void followTrack(GbtsEdge& pS, GbtsEdgeState& output);

 protected:
  /// Propagate edge state
  /// @param pS Edge to propagate from
  /// @param ts Edge state to update
  void propagate(GbtsEdge& pS, GbtsEdgeState& ts);

  /// Update edge state with edge
  /// @param pS Edge to update with
  /// @param ts Edge state to update
  /// @return Success flag
  bool update(GbtsEdge& pS, GbtsEdgeState& ts);

  /// Get layer type from layer index
  /// @param l Layer index
  /// @return Layer type
  std::uint32_t getLayerType(std::uint32_t l);

  /// Geometry layer description
  const std::vector<TrigInDetSiLayer>& m_geo;

  /// Edge storage
  std::vector<GbtsEdge>& m_segStore;

  /// State vector
  std::vector<GbtsEdgeState*> m_stateVec;

  /// State storage array
  std::array<GbtsEdgeState, GbtsMaxEdgeState> m_stateStore;

  /// Global state counter
  std::uint32_t m_globalStateCounter{0};

  /// Configuration for seed finder
  const SeedFinderGbtsConfig& m_config;
};

}  // namespace Acts::Experimental
