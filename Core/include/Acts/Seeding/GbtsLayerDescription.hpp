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

/// Identifier of a GBTS layer as the experiment numbers it, which is what the
/// layer connection table is written in terms of. Sparse, and an index into
/// nothing.
using GbtsLayerId = std::uint32_t;

/// Position of a layer in the GBTS geometry, assigned in the order the layer
/// descriptions are handed to it. This is what indexes the geometry and what a
/// node carries.
using GbtsLayerIndex = std::uint32_t;

/// Where a GBTS layer sits, which fixes its coordinate convention.
enum class GbtsLayerType : std::uint8_t { Barrel = 0, Endcap = 1 };

/// Sensor technology of a GBTS layer, independent of GbtsLayerType.
enum class GbtsLayerTechnology : std::uint8_t { Pixel = 0, Strip = 1 };

/// Lightweight layer description for GBTS geometry.
struct GbtsLayerDescription final {
  /// Combined subdetector ID.
  GbtsLayerId id{};
  /// Layer type (barrel or endcap).
  GbtsLayerType type{};
  /// Sensor technology of the layer.
  GbtsLayerTechnology technology{GbtsLayerTechnology::Pixel};
  /// Reference coordinate (r for barrel, z for endcap).
  float refCoord{};
  /// Minimum boundary coordinate.
  float minBound{};
  /// Maximum boundary coordinate.
  float maxBound{};
};

}  // namespace Acts::Experimental
