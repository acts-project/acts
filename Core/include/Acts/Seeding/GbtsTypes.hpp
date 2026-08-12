// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"

#include <array>
#include <cstddef>
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

/// Number of values per row of the tau lookup table: the cluster width the row
/// is keyed on, followed by two pairs of tau bounds.
static constexpr std::size_t kGbtsTauLookupTableColumns = 5;

/// Per-cluster-width tau bounds, one row per 0.05 mm of cluster width. The
/// bounds are trained offline; the seeder only reads the table back.
using GbtsTauLookupTable =
    std::vector<std::array<float, kGbtsTauLookupTableColumns>>;

}  // namespace Acts::Experimental
