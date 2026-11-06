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

  /// Transverse state in the frame of the head edge: y', dy'/dx', curvature
  std::array<float, 3> x{};
  /// Longitudinal state: z, dz/ds
  std::array<float, 2> y{};
  /// Transverse covariance, upper triangle by rows: 00, 01, 02, 11, 12, 22
  std::array<float, 6> cx{};
  /// Longitudinal covariance, upper triangle by rows: 00, 01, 11
  std::array<float, 3> cy{};
  /// Position of the state: x' and r
  float refX{};
  float refY{};
  /// Cosine and sine of the direction of the head edge, defining x'
  float c{};
  float s{};
};

}  // namespace Acts::Experimental::detail
