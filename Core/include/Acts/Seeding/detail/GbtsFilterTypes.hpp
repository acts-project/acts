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
  /// @param varianceX Initial variance of the transverse state
  /// @param varianceY Initial variance of the longitudinal state
  void initialize(const GbtsEdge& pS, const GbtsNodeView& nodeView,
                  const std::array<float, 3>& varianceX,
                  const std::array<float, 2>& varianceY);

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

/// State of the search for the best fitted chain from an edge.
struct GbtsFilterState final {
  /// Minimum number of edges of a complete chain.
  std::uint32_t minChainLength{0};

  std::vector<std::uint32_t> chain;

  std::vector<std::uint32_t> bestChain;

  /// Seed quality of the best chain, lower is better.
  float bestQuality{0};

  GbtsEdgeState bestState;
};

}  // namespace Acts::Experimental::detail
