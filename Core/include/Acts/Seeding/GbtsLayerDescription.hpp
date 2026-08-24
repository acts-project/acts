// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstdint>

namespace Acts::Experimental {

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
  /// Half width in reference direction
  float halfRefWidth = {};
  /// Minimum boundary coordinate.
  float minBound{};
  /// Maximum boundary coordinate.
  float maxBound{};
};

}  // namespace Acts::Experimental
