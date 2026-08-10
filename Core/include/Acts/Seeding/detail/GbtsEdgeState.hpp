// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/GbtsGraphTypes.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace Acts::Experimental::detail {

/// Per-edge tracking state used by the GBTS filter.
struct GbtsEdgeState final {
 public:
  GbtsEdgeState() = default;

  /// Constructor with initialization flag
  /// @param f Initialization flag
  explicit GbtsEdgeState(bool f) : initialized(f) {}

  /// Initialize from edge
  /// @param pS Edge to initialize from
  /// @param nodeView View of the node positions and layers
  void initialize(const GbtsEdge& pS, const GbtsNodeView& nodeView);

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

/// Maximum number of edge states
static constexpr std::uint32_t kGbtsMaxEdgeStates = 2500;

/// State for the tracking filter, containing edge states and a global counter.
struct GbtsFilterState final {
  /// State vector
  std::vector<GbtsEdgeState*> stateVec;

  /// State storage array
  std::array<GbtsEdgeState, kGbtsMaxEdgeStates> stateStore{};

  /// Global state counter
  std::uint32_t globalStateCounter{0};
};

}  // namespace Acts::Experimental::detail
