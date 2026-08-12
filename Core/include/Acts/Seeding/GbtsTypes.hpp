// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"

#include <cstdint>
#include <vector>

namespace Acts::Experimental {

/// Index of a graph node inside GbtsNodeStorage. Nodes are stored ordered by
/// (eta bin, phi), so all nodes of an eta bin form a contiguous range.
using GbtsNodeIndex = SpacePointIndex;

/// Sentinel for an unset graph node index.
static constexpr GbtsNodeIndex kGbtsNodeIndexInvalid = kSpacePointIndexInvalid;

/// GBTS layer types
enum class GbtsLayerType { Barrel = 0, Endcap = 1 };

/// Lightweight layer description for GBTS geometry.
struct GbtsLayerDescription final {
  /// Combined subdetector ID.
  std::int32_t id{};
  /// Layer type (barrel or endcap).
  GbtsLayerType type{};
  /// Reference coordinate (z for barrel, r for endcap).
  float refCoord{};
  /// Minimum boundary coordinate.
  float minBound{};
  /// Maximum boundary coordinate.
  float maxBound{};
};

/// Accepted |cot(theta)| range for one cluster width bin. A cluster near the
/// module edge may be shortened, hence the second pair.
struct GbtsTauBounds final {
  /// Minimum accepted |cot(theta)|.
  float minTau{};
  /// Maximum accepted |cot(theta)|. Negative means undertrained: do not cut.
  float maxTau{};
  /// Minimum accepted |cot(theta)| near the module edge.
  float minTauNearEdge{};
  /// Maximum accepted |cot(theta)| near the module edge. Negative as above.
  float maxTauNearEdge{};
};

/// Tau bounds per cluster width, one entry per 0.05 mm, indexed not searched.
using GbtsTauLookupTable = std::vector<GbtsTauBounds>;

}  // namespace Acts::Experimental
