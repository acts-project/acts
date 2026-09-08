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

/// Layer id as the experiment numbers it (80000, 81000, ...). Sparse, and an
/// index into nothing.
using GbtsExperimentLayerId = std::uint32_t;

/// Position of a layer within one GbtsGeometry, dense from zero. Ask
/// GbtsGeometry::layerIndex for it rather than deriving it.
using GbtsLayerIndex = std::uint16_t;

/// Where a GBTS layer sits, which fixes its coordinate convention.
enum class GbtsLayerType : std::uint8_t { Barrel = 0, Endcap = 1 };

/// Sensor technology of a GBTS layer, independent of GbtsLayerType.
enum class GbtsLayerTechnology : std::uint8_t { Pixel = 0, Strip = 1 };

/// Lightweight layer description for GBTS geometry.
struct GbtsLayerDescription final {
  /// Combined subdetector ID.
  GbtsExperimentLayerId id{};
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
  /// Position in the inside-out ordering of the pixel barrel layers.
  std::int32_t barrelOrder{-1};
  /// Cut the layer's nodes against the z0 range.
  /// By default, only the innermost barrel layer.
  bool cutOnZ0Range{false};
  /// Apply `matchBeforeCreate` to the layer.
  /// By default, the first two barrel layers.
  bool matchBeforeCreate{false};
};

}  // namespace Acts::Experimental
