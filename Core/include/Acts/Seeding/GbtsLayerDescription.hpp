// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <cstdint>
#include <vector>

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
  /// Thickness in the reference direction
  /// max r - min r for the barrel, max z - min z for endcap
  float layerThickness{};
  /// Minimum boundary coordinate.
  float minBound{};
  /// Maximum boundary coordinate.
  float maxBound{};
  /// Position of the layer in the inside-out ordering of the pixel barrel,
  /// -1 for every other layer. Derived from `refCoord` when left unset; set it
  /// on every pixel barrel layer or on none of them.
  std::int32_t barrelOrder{-1};
};

/// A GBTS layer and the detector surfaces it is made of
struct GbtsLayerConfig final {
  /// Id of the layer, as the layer connections refer to it
  GbtsExperimentLayerId id{};
  /// Whether the layer is a barrel or an endcap layer
  GbtsLayerType type{};
  /// Sensor technology of the layer
  GbtsLayerTechnology technology{};
  /// Surfaces of the layer
  std::vector<GeometryIdentifier> surfaces;
};

}  // namespace Acts::Experimental
