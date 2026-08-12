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

  bool initialized{false};

  /// Score for comparison
  float j{};

  std::vector<GbtsEdge*> vs;

  std::array<float, 3> x{};
  std::array<float, 2> y{};
  std::array<std::array<float, 3>, 3> cx{};
  std::array<std::array<float, 2>, 2> cy{};
  float refX{};
  float refY{};
  float c{};
  float s{};
};

/// Maximum number of edge states
static constexpr std::uint32_t kGbtsMaxEdgeStates = 2500;

/// State for the tracking filter, containing edge states and a global counter.
struct GbtsFilterState final {
  std::vector<GbtsEdgeState*> stateVec;

  std::array<GbtsEdgeState, kGbtsMaxEdgeStates> stateStore{};

  std::uint32_t globalStateCounter{0};
};

}  // namespace Acts::Experimental::detail
