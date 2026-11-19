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

namespace Acts::Experimental::detail {

/// Per-edge tracking state used by the GBTS filter.
struct GbtsEdgeState final {
 public:
  GbtsEdgeState() = default;

  /// Initialize from edge
  /// @param pS Edge to initialize from
  /// @param nodeView View of the node positions and layers
  /// @param varianceX Initial variance of the transverse state
  /// @param varianceY Initial variance of the longitudinal state
  void initialize(const GbtsEdge& pS, const GbtsNodeView& nodeView,
                  const std::array<float, 3>& varianceX,
                  const std::array<float, 2>& varianceY);

  /// Score for comparison
  float j{};

  std::array<float, 3> x{};
  std::array<float, 2> y{};
  std::array<std::array<float, 3>, 3> cx{};
  std::array<std::array<float, 2>, 2> cy{};
  float refX{};
  float refY{};
  float c{};
  float s{};
};

}  // namespace Acts::Experimental::detail
